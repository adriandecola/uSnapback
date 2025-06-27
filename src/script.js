/*
File:           script.js
Description:    Main JavaScript file for interactivity, functionality, and processing.
Author:         Adrian deCola
Relative Path:  uSnapback/src/script.js
*/

/*****************************************************************************************/
/*************************************** Constants ***************************************/
/*****************************************************************************************/
const NUCLEOTIDE_COMPLEMENT = { A: 'T', T: 'A', C: 'G', G: 'C' };
const STRONG_NUCLEOTIDE_MISMATCH = { A: 'G', G: 'A', C: 'C', T: 'T' };
const VALID_BASES = new Set(['A', 'T', 'C', 'G']);
const SNV_BASE_BUFFER = 4; // The number of matched bases required on either end of a mismatched SNV
const INNER_LOOP_NUMBER_OF_STRONG_BASE_MISMATCHES_REQUIRED = 2;
const END_OF_STEM_NUMBER_OF_STRONG_BASE_MISMATCHES_REQUIRED = 2;
const MINIMUM_TARGET_SNAPBACK_MELTING_TEMP = 40;
const MAXIMUM_TARGET_SNAPBACK_MELTING_TEMP = 80;
const MIN_LOOP_LEN = 6;
const MIN_PRIMER_LEN = 12;
const TM_DECIMAL_PLACES = 2;
const API_TOKEN = 1;
// Chemisty parameters ************** name these better ************
const MG = 2.2;
const MONO = 20.0;
const T_PARAM = 'UnifiedSantaLucia';
const SALT_CALC_TYPE = 'bpdenominator';
const O_TYPE = 'oligo';
const CONC = 0.5;
const LIMITING_CONC = 0.5;

/*****************************************************************************************/
/************************************ Primary Function ***********************************/
/*****************************************************************************************/
/**
 * Creates a snapback primer by:
 *  1. Evaluating both strands to decide which primer receives the tail and
 *     whether the tail pairs to the wild or variant allele (choose the option
 *     with the largest wild/variant Tm difference in an initial small stem).
 * 	   This initial stem is built by adding `SNV_BASE_BUFFER` bases on each end
 *     of the single nucleotide variant (SNV).
 *  2. Extending that initial stem outward, respecting primer annealing locations,
 *     until the wild-type stem Tm approaches `targetSnapMeltTemp`, or no
 *     further growth is possible.  This temperature must exceed `minSnapbackMeltTemp`.
 *  3. Building the final 5'→3' sequence:
 *        snapback-tail • inner-loop mismatches • stem (with chosen SNV base) • primer.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Assumptions
 * ──────────────────────────────────────────────────────────────────────────
 * - targetSeqStrand is a valid uppercase DNA string given 5'→3'.
 * - primerLen and compPrimerLen ≥ {MIN_PRIMER_LEN}.
 * - The SNV is ≥ {SNV_BASE_BUFFER} bases away from both primers.
 * - Stem Tm values are computed with the Santa Lucia nearest-neighbour model
 *   via dna-utah.org.
 * - We want to miminize the loop length to the primer on which the snapback tail is on
 *   plus the two strong mismatched in the inner loop so it doesn't zip
 * - We want to keep the SNV in the middle of the stem or as close to centered as we
 *   can, as we build the snapback stem
 * - Extension can occur on the snapback's complement, on the side of the stem does not contain
 *   the loop, if the polymerase can easily attach to that location. A 2 base pair, strong mismatch,
 *   is added after the 5' end of the stem (on the snapback primer) help avoid this extension
 *   on its complement snapback
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Possible improvements
 * ──────────────────────────────────────────────────────────────────────────
 * - I could represent DNA sequences as an array of 2 bit encoded neucleotides. This could save some memory,
 *   albeit minimal. It would be a fun thing to code
 * 		- Then DNA functions could be implemented as methods
 * - I could implement the SantaLucia melting temperature calculations for the STEM in JavaScript in this file
 *   as it's own function
 * - I could add in some temperature correction for the dangling end/end mismatch for the snapback primer on
 *   on the end of the stem (the strong mismatches that prevent extension on the complementary primer)
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Type definitions
 * ──────────────────────────────────────────────────────────────────────────
 * @typedef {Object} SNVSite
 * @property {number} 				index				0-based position of the SNV on `targetSeqStrand`
 * @property {string}				variantBase			Variant base ('A', 'C', 'G', or 'T')
 *
 * @typedef {Object} SnapbackMeltingTm
 * @property {number} 				wildTm     			Calculated Tm (°C) of the snapback on the wilt type allele
 * @property {number} 				variantTm  			Calculated Tm (°C) of the snapback on the variant type allele
 *
 * @typedef {Object} SnapbackPrimerResult
 * @property {string}				snapbackSeq			Entire snapback primer written 5' → 3'
 *                                                   	(tail → primer).
 * @property {boolean}				isOnTargetPrimer	true if tail is appended to the forward primer, i.e. the
 * 														primer on `targetSeqStrand`; false if it is appended
 * 														to the reverse primer.
 * @property {boolean}				matchesWild			true if the snapback base at the SNV matches
 *                                       				the wild-type allele
 * @property {SnapbackMeltingTm}	snapbackMeltingTm	Object holding wild/variant snapback Tm values.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Parameters, Returns, and Errors
 * ──────────────────────────────────────────────────────────────────────────
 * @param {string}				targetSeqStrand					A strand of the full DNA sequence to design the snapback primer for.
 * @param {number}				primerLen						The length of the forward primer, primer that shares the nucleotides with the
 * 																strand denoted by the `targetSeqStrand`
 * @param {number}				compPrimerLen					The length of the complementary primer on the complementary strand to
 * 																`targetSeqStrand`
 * @param {SNVSite}				snvSite							An object representing the single nucleotide variant site
 * @param {number}				targetSnapMeltTemp	The desired snapback melting temperature for the wild type allele
 *
 * @returns {Promise<SnapbackPrimerResult>}						An object representing the formed snapback primer and is specification
 *
 * @throws {Error} 												If any input is invalid, the SNV is too close to a primer,
 *             													or an acceptable stem cannot be constructed.
 */
async function createSnapback(
	targetSeqStrand,
	primerLen,
	compPrimerLen,
	snvSite,
	targetSnapMeltTemp
) {
	//──────────────────────────────────────────────────────────────────────────//
	//							Parameter Checking								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1. targetSeqStrand
	if (!isValidDNASequence(targetSeqStrand)) {
		throw new Error(
			`Invalid DNA sequence: ${targetSeqStrand}. ` +
				`Must be a non-empty uppercase string containing only characters A, T, C, and/or G..`
		);
	}

	// Validate targetSnapMeltTemp is a positive number
	if (
		typeof targetSnapMeltTemp !== 'number' ||
		!Number.isFinite(targetSnapMeltTemp) ||
		targetSnapMeltTemp < 0
	) {
		throw new Error(`targetSnapMeltTemp must be a positive, finite number`);
	}

	// Validate primerLen and compPrimerLen are positive integers
	for (const [name, val] of [
		['primerLen', primerLen],
		['compPrimerLen', compPrimerLen],
	]) {
		if (typeof val !== 'number' || !Number.isFinite(val)) {
			throw new Error(`${name} must be a finite number`);
		}
		if (val < 0) {
			throw new Error(`${name} must be non-negative`);
		}
		if (!Number.isInteger(val)) {
			throw new Error(`${name} must be an integer`);
		}
		if (val < MIN_PRIMER_LEN) {
			throw new Error(`${name} must be at lease ${MIN_PRIMER_LEN}`);
		}
	}

	// Validate snvSite structure and values
	if (
		typeof snvSite !== 'object' ||
		snvSite == null ||
		Array.isArray(snvSite)
	) {
		throw new Error(`snvSite must be a non-null object`);
	}
	const expectedKeys = new Set(['index', 'variantBase']);
	for (const key of Object.keys(snvSite)) {
		if (!expectedKeys.has(key)) {
			throw new Error(`Unexpected key in snvSite: "${key}"`);
		}
	}
	if (!('index' in snvSite) || !('variantBase' in snvSite)) {
		throw new Error(`snvSite must contain both "index" and "variantBase"`);
	}
	if (
		typeof snvSite.index !== 'number' ||
		!Number.isInteger(snvSite.index) ||
		snvSite.index < 0 ||
		snvSite.index >= targetSeqStrand.length
	) {
		throw new Error(
			`snvSite.index must be an integer from 0 to ${
				targetSeqStrand.length - 1
			}`
		);
	}
	if (
		typeof snvSite.variantBase !== 'string' ||
		snvSite.variantBase.length !== 1 ||
		!VALID_BASES.has(snvSite.variantBase)
	) {
		throw new Error(
			`snvSite.variantBase must be a single character: A, T, C, or G`
		);
	}
	//──────────────────────────────────────────────────────────────────────────//
	//								Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1) Make sure the SNV is not too close to either end of the primer that a proper stem cannot be formed
	if (
		snvTooCloseToPrimer(
			snvSite.index,
			primerLen,
			compPrimerLen,
			targetSeqStrand.length
		)
	) {
		throw new Error(
			`SNV at index ${snvSite.index} is too close to one of the primers. ` +
				`There must be at least 3 bases on each side of the SNV that do not ` +
				`over lap with either primer location to form a valid stem.`
		);
	}

	// 2) Calculate the initial melting temperature differences for variants and choose the primer and base match (variant or wild)
	//	  with the largest differnce as snapback
	// 	  Note: It is assumed that the melting temperature difference will scale(decrease) stem length is increased, but that the
	//	  snapback with the greatest temperature difference will remain the same
	const { useTargetStrand, snapbackBaseAtSNV, matchesWild } =
		await useTargetStrandsPrimerForComplement(targetSeqStrand, snvSite);

	// 3) Assigning variables in terms of the primer strand to use as the snapback
	var targetStrandSeqSnapPrimerRefPoint;
	var snvSiteSnapPrimerRefPoint;
	if (useTargetStrand) {
		targetStrandSeqSnapPrimerRefPoint = targetSeqStrand;
		snvSiteSnapPrimerRefPoint = snvSite;
	} else {
		targetStrandSeqSnapPrimerRefPoint = reverseComplement(targetSeqStrand);
		snvSiteSnapPrimerRefPoint = revCompSNV(snvSite);
	}

	// 4) Calculating the stem
	const { stemLoc, meltingTemps } = createStem(
		targetStrandSeqSnapPrimerRefPoint,
		snvSiteSnapPrimerRefPoint,
		primerLensSnapPrimerRefPoint,
		snapbackBaseAtSNV,
		matchesWild,
		targetSnapMeltTemp
	);

	// 5) Create a final snapback.
	const snapback = buildFinalSnapback(
		targetStrandSeqSnapPrimerRefPoint,
		snvSiteSnapPrimerRefPoint,
		primerLenstemLoc,
		snapbackBaseAtSNV
	);
}

/*****************************************************************************************/
/********************************** Secondary Functions **********************************/
/*****************************************************************************************/

/**
 * Decides which strand (target vs. complementary) should get the snapback primer,
 * and whether the snapback should match the wild-type base or the single variant base at the SNV.
 *
 * @typedef {Object} snapbackDecision
 *    @property {boolean} useTargetStrand - True if the snapback tail should be added to the target strand.
 *    @property {string} snapbackBaseAtSNV - The base (e.g., "A", "G") used in the snapback at the SNV site.
 *    @property {boolean} matchesWild - True if the snapback matches the wild-type base.
 *
 * @param {string} targetSeqStrand     - Full target strand (5'→3')
 * @param {Object} snvSite             - { index: number, variantBase: string }
 *
 * @returns {Promise<{ snapbackDecision}>} - The decision object for which primer should be the snapback and what it should match
 * 											 (wild or variant type)
 */
async function useTargetStrandsPrimerForComplement(targetSeqStrand, snvSite) {
	//──────────────────────────────────────────────────────────────────────────//
	//							Parameter Checking								//
	//──────────────────────────────────────────────────────────────────────────//
	// 1. Validate targetSeqStrand
	if (!isValidDNASequence(targetSeqStrand)) {
		throw new Error(`Invalid targetSeqStrand: "${targetSeqStrand}"`);
	}

	// 2. Validate SNV site object
	// Checks that SNV sites are JavaScript objects
	if (
		typeof snvSite !== 'object' ||
		snvSite == null ||
		Array.isArray(snvSite)
	) {
		throw new Error(
			`snvSite must be an object with 'index' and 'variantBase' fields.`
		);
	}
	// Checks proper keys, and only proper keys are part of SNV site Javascript objects
	const expectedKeys = new Set(['index', 'variantBase']);
	const actualKeys = new Set(Object.keys(snvSite));
	for (const key of actualKeys) {
		if (!expectedKeys.has(key)) {
			throw new Error(`snvSite contains unexpected key: "${key}"`);
		}
	}
	if (!('index' in snvSite) || !('variantBase' in snvSite)) {
		throw new Error(
			`snvSite must include both "index" and "variantBase" keys.`
		);
	}
	// Validate values of keys
	if (
		!Number.isInteger(snvSite.index) ||
		snvSite.index < 0 ||
		snvSite.index >= targetSeqStrand.length
	) {
		throw new Error(
			`snvSite.index must be a valid integer from 0 to ${
				targetSeqStrand.length - 1
			}`
		);
	}
	if (
		typeof snvSite.variantBase !== 'string' ||
		snvSite.variantBase.length !== 1 ||
		!VALID_BASES.has(snvSite.variantBase)
	) {
		throw new Error(
			`snvSite.variantBase must be a single character: A, T, C, or G`
		);
	}

	// 3) Ensure the sequence is long enough around the SNV site to accommodate the initial stem
	if (
		snvSite.index < SNV_BASE_BUFFER ||
		snvSite.index > targetSeqStrand.length - SNV_BASE_BUFFER - 1
	) {
		throw new Error(
			`Cannot create initial stem: SNV at index ${snvSite.index} is too close to the ends of the sequence. ` +
				`Need ${SNV_BASE_BUFFER} matched bases on both sides. Sequence length is ${targetSeqStrand.length}.`
		);
	}

	//──────────────────────────────────────────────────────────────────────────//
	//								Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1) Build reverse complement of the target sequence's strand
	const revCompTargetSeqStrand = reverseComplement(targetSeqStrand);

	// 2) Build the SNV Site object for the reveres complement strand
	const revCompSnvSite = revCompSNV(snvSite, targetSeqStrand.length);

	// 3) Adds {SNV base buffer} matching neucleotides on each end of SNV to create the initial stem
	const initStemLoc = {
		start: snvSite.index - SNV_BASE_BUFFER,
		end: snvSite.index + SNV_BASE_BUFFER,
	};
	const compInitStemLoc = {
		start: revCompSnvSite.index - SNV_BASE_BUFFER,
		end: revCompSnvSite.index + SNV_BASE_BUFFER,
	};

	// 4) Slice out the "init stem" region from each strand
	const targetInitStem = targetSeqStrand.slice(
		initStemLoc.start,
		initStemLoc.end + 1
	);
	const compInitStem = revCompTargetSeqStrand.slice(
		compInitStemLoc.start,
		compInitStemLoc.end + 1
	);

	// 5) Identify the wild-type base on each strand
	//    Should be complementary
	const targetWildBase = targetSeqStrand[snvSite.index];
	const compWildBase = revCompTargetSeqStrand[revCompSnvSite.index];

	// 6) SNV is at position {SNV_BASE_BUFFER} in these 2*{SNV_BASE_BUFFER}+1 slices
	const mismatchPos = SNV_BASE_BUFFER;

	// 7) Evaluate Tm differences for snapback tail on target strand
	const targetScenario = await evaluateSnapbackMatchingOptions(
		targetInitStem,
		mismatchPos,
		targetWildBase,
		snvSite.variantBase
	);

	// 8) Evaluate Tm differences for snapback tail on complementary strand
	const compScenario = await evaluateSnapbackMatchingOptions(
		compInitStem,
		mismatchPos,
		compWildBase,
		revCompSnvSite.variantBase
	);

	// 9) Compare which scenario yields the bigger Tm difference
	if (targetScenario.bestDifference > compScenario.bestDifference) {
		return {
			useTargetStrand: true,
			snapbackBaseAtSNV: targetScenario.bestSnapbackBase,
			matchesWild: targetScenario.matchesWild,
		};
	} else {
		return {
			useTargetStrand: false,
			snapbackBaseAtSNV: compScenario.bestSnapbackBase,
			matchesWild: compScenario.matchesWild,
		};
	}
}

/**
 * Evaluates the potential snapback Tm differences for a given stem slice
 * when there's only a single variant base.
 *
 * We consider two scenarios:
 *  1) Snapback matches the wild base => the variant is the mismatch.
 *  2) Snapback matches the variant base => the wild base is the mismatch.
 *
 * Returns whichever scenario yields the largest absolute difference from the
 * “matched wild” Tm.
 *
 * @param {string} initStem     - The sub-sequence {2*SNV_BASE_BUFFER+1} nucleotides long, containing the SNV at position `mismatchPos`
 * 								  assumed to be passed in with the wild base.
 * @param {number} mismatchPos  - Usually {SNV_BASE_BUFFER} (SNV in center)
 * @param {string} wildBase     - The wild-type base at that SNV
 * @param {string} variantBase  - The single variant base
 *
 * @returns {Promise<{ bestSnapbackBase: string, bestDifference: number, matchesWild: boolean }>}
 *   bestSnapbackBase => the base to use in the snapback (either complementary to the wild or variant base)
 *   bestDifference => the numeric Tm difference from the wild scenario
 *   matchesWild => true if the snapback tail matches the wild type at the SNV location
 */
async function evaluateSnapbackMatchingOptions(
	initStem,
	mismatchPos,
	wildBase,
	variantBase
) {
	//──────────────────────────────────────────────────────────────────────────//
	//							Parameter Checking								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1. Validate initStem
	if (!isValidDNASequence(initStem)) {
		throw new Error(`Invalid initStem sequence: "${initStem}"`);
	}

	// 2. Validate mismatchPos
	if (
		typeof mismatchPos !== 'number' ||
		!Number.isInteger(mismatchPos) ||
		mismatchPos < 0 ||
		mismatchPos >= initStem.length
	) {
		throw new Error(
			`mismatchPos must be an integer in range [0, ${
				initStem.length - 1
			}]. Received: ${mismatchPos}`
		);
	}

	// 3. Validate wildBase
	if (
		typeof wildBase !== 'string' ||
		wildBase.length !== 1 ||
		!VALID_BASES.has(wildBase)
	) {
		throw new Error(
			`wildBase must be a single character from "A", "T", "C", or "G". Received: "${wildBase}"`
		);
	}
	// Check that initStem[mismatchPos] === wildBase
	if (initStem[mismatchPos] !== wildBase) {
		throw new Error(
			`Mismatch position ${mismatchPos} in initStem does not contain wildBase. ` +
				`Found "${initStem[mismatchPos]}", expected "${wildBase}".`
		);
	}

	// 4. Validate variantBase
	if (
		typeof variantBase !== 'string' ||
		variantBase.length !== 1 ||
		!VALID_BASES.has(variantBase)
	) {
		throw new Error(
			`variantBase must be a single character from "A", "T", "C", or "G". Received: "${variantBase}"`
		);
	}

	// 5. Ensure variantBase differs from wildBase
	if (variantBase === wildBase) {
		throw new Error(
			`variantBase and wildBase must differ to represent a true SNV. Both were "${wildBase}".`
		);
	}

	//──────────────────────────────────────────────────────────────────────────//
	//								Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1) Tm with wild base, matched
	const wildMatchTm = await getStemTm(initStem);

	// 2) Tm with variant base, matched
	// Creating initial stem for the variant sequence
	const variantInitStem =
		initStem.slice(0, mismatchPos) +
		variantBase +
		initStem.slice(mismatchPos + 1);

	const variantMatchTm = await getStemTm(variantInitStem);

	// 2) Scenario A: Wild-matching snapback (mis)matched to the variant target sequence
	const wildSnapbacktoVariantMismatchObj = {
		position: mismatchPos,
		// At the mismatch, the snapback nucleotide will be the complement of the wild type base
		type: NUCLEOTIDE_COMPLEMENT[wildBase],
	};
	const wildSnapbacktoVariantTm = await getStemTm(
		variantInitStem,
		wildSnapbacktoVariantMismatchObj
	);
	const wildSnapbackTmDiff = Math.abs(wildMatchTm - wildSnapbacktoVariantTm);

	// 3) Scenario B: Variant-matching snapback (mis)matches to the wild target sequence
	const variantSnapbacktoWildMismatchObj = {
		position: mismatchPos,
		// At the mismatch, the snapback nucleotide will be the complement of the variant type base
		type: NUCLEOTIDE_COMPLEMENT[variantBase],
	};
	const variantSnapbacktoWildTm = await getStemTm(
		initStem,
		variantSnapbacktoWildMismatchObj
	);
	const variantSnapbackTmDiff = Math.abs(
		variantMatchTm - variantSnapbacktoWildTm
	);

	// 4) Pick whichever scenario yields the larger difference
	if (wildSnapbackTmDiff > variantSnapbackTmDiff) {
		// Scenario A wins
		// Snapback should match wild
		return {
			bestSnapbackBase: NUCLEOTIDE_COMPLEMENT[wildBase],
			bestDifference: wildSnapbackTmDiff,
			matchesWild: true, // <-- ADDED
		};
	} else {
		// Scenario B wins
		// Snapback should match variant
		return {
			bestSnapbackBase: NUCLEOTIDE_COMPLEMENT[variantBase],
			bestDifference: variantSnapbackTmDiff,
			matchesWild: false,
		};
	}
}

/**
 * Calculates the melting temperature of a snapback stem via server API (dna-utah.org).
 *If the mismatch object is passed in, it calculates the mismatch Tm
 *
 * @typedef {Object} Mismatch
 * @property {number} position - Index in the sequence where the mismatch occurs.
 * @property {string} type     - The base on the opposite strand (e.g. "A") for the mismatch.
 *
 * @param {string}  seq             - The reference (matched) sequence (5'→3').
 * @param {Mismatch} [mismatch]     - Optional mismatch specification.
 * @returns {Promise<number>}       - The Tm value in °C.
 * @throws {Error}                  - If inputs are invalid, or the response can’t be parsed.
 */
async function getStemTm(seq, mismatch) {
	//──────────────────────────────────────────────────────────────────────────//
	//							Parameter Checking								//
	//──────────────────────────────────────────────────────────────────────────//
	if (!isValidDNASequence(seq)) {
		throw new Error(`Invalid, empty, or non-string DNA sequence: "${seq}"`);
	}
	if (mismatch !== undefined && mismatch !== null) {
		// it could be passed in as null and we'll ignore it
		if (typeof mismatch !== 'object' || Array.isArray(mismatch)) {
			throw new Error(
				`Mismatch must be an object. Received: ${typeof mismatch}`
			);
		}
		// Check for extra keys in mismatch object
		const allowedMismatchKeys = new Set(['position', 'type']);
		for (const key of Object.keys(mismatch)) {
			if (!allowedMismatchKeys.has(key)) {
				throw new Error(
					`Mismatch object contains unexpected key: "${key}"`
				);
			}
		}
		if (!('position' in mismatch) || !('type' in mismatch)) {
			throw new Error(
				`Mismatch object missing required keys "position" and/or "type". Received: ${JSON.stringify(
					mismatch
				)}`
			);
		}
		if (
			typeof mismatch.position !== 'number' ||
			!Number.isInteger(mismatch.position) ||
			mismatch.position < 0 ||
			mismatch.position >= seq.length
		) {
			throw new Error(
				`Mismatch position ${mismatch.position} is invalid or out of bounds for sequence of length ${seq.length}`
			);
		}
		if (
			typeof mismatch.type !== 'string' ||
			mismatch.type.length !== 1 ||
			!VALID_BASES.has(mismatch.type)
		) {
			throw new Error(
				`Mismatch type "${mismatch.type}" must be one of "A", "T", "C", "G"`
			);
		}
	}

	//──────────────────────────────────────────────────────────────────────────//
	//								Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//

	// 2) Build the mismatch sequence the way the API endpoint wants it, if its provided
	let mismatchSeq = null;
	if (mismatch) {
		try {
			mismatchSeq = buildMismatchSequenceForAPI(seq, mismatch);
		} catch (err) {
			// Re-throw error
			throw new Error(
				`Failed to build mismatch sequence: ${err.message}`
			);
		}
	}

	// 3) Build query URL
	let baseUrl = 'https://dna-utah.org/tm/snaprequest.php';
	baseUrl += `?mg=${MG}`;
	baseUrl += `&mono=${MONO}`;
	baseUrl += `&seq=${seq.toLowerCase()}`;
	baseUrl += `&tparam=${T_PARAM}`;
	baseUrl += `&saltcalctype=${SALT_CALC_TYPE}`;
	baseUrl += `&otype=${O_TYPE}`;
	baseUrl += `&concentration=${CONC}`;
	baseUrl += `&limitingconc=${LIMITING_CONC}`;
	baseUrl += `&decimalplaces=${TM_DECIMAL_PLACES}`;
	baseUrl += `&token=${API_TOKEN}`;
	// If a mismatch is passed in this will add the correct mismatch sequence
	if (mismatch) {
		baseUrl += `&mmseq=${mismatchSeq}`;
	}

	// 5) Fetch and parse response
	const response = await fetch(baseUrl);
	if (!response.ok) {
		throw new Error(
			`Network error: ${response.status} - ${response.statusText}`
		);
	}
	// Parses the body of the response as text
	const rawHtml = await response.text();

	// Parsing the correct tm value (tm or mmtm)
	let tmValue;
	if (!mismatch) {
		tmValue = parseTmFromResponse(rawHtml);
	} else {
		tmValue = parseTmFromResponse(rawHtml, true);
	}

	// Throw error if tmValue is not found
	if (tmValue === null) {
		throw new Error(
			'No <tm> element found or invalid numeric value in server response.'
		);
	}

	// 6) Return the numeric Tm
	return tmValue;
}

/**
 * Constructs the snapback stem on the chosen strand/primer so that
 * the melting temperature of the wild type allele gets as close as possible
 * to the desired snapback melting temperature. It does this while keeping the loop length
 * mimimized to the primer length plus the number of strong inner loop mismatches required.
 *
 * @typedef {Object} SNVSiteRefPoint
 * @property {number} index - The 0-based index of the SNV on this strand's coordinate system.
 * @property {string} variantBase - The variant base at this position (must be "A", "T", "C", or "G").
 *
 * @typedef {Object} PrimerLensRefPoint
 * @property {number} primerLen - The length of the primer on this strand (snapback primers complementary strand)
 * @property {number} compPrimerLen - The length of the complementary primer on the opposite strand
 *
 * @typedef {Object} MeltingTemp
 * @property {number} wildTm - The snapback's melting temperature for the wild type allele
 * @property {number} variantTm - The snapback's melting temperature for the variant type allele
 *
 * @typedef {Object} StemLoc
 * @property {number} start - The start index (inclusive) of the stem in the strand (0-based).
 * @property {number} end - The end index (inclusive) of the stem in the strand (0-based).
 *
 * @typedef {Object} CreateStemReturn
 * @property {StemLoc} stemLoc - The finalized stem location for the snapback.
 * @property {MeltingTemp} meltingTemps - Calculated Tm values for both wild-type and variant.
 * @property {string} snapbackBaseAtSNV - The base (complement of wild or variant) to use in the snapback at SNV site.
 *
 * @param {string} targetStrandSeqSnapPrimerRefPoint - DNA sequence (5'→3') of the strand used for the snapback.
 * @param {SNVSiteRefPoint} snvSiteSnapPrimerRefPoint - SNV information referenced in this strand's frame.
 * @param {PrimerLensRefPoint} primerLensSnapPrimerRefPoint - Object containing both the primer and complementary primer lengths.
 * @param {string} snapbackBaseAtSNV - The base to place in the snapback tail at the SNV position.
 * @param {boolean} matchesWild - Whether the snapback tail matches the wild-type allele.
 * @param {number} targetSnapMeltTemp - Desired melting temperature (°C) for the wild-type stem.
 *
 * @returns {CreateStemReturn} Finalized stem location and melting temperature information.
 */
async function createStem(
	targetStrandSeqSnapPrimerRefPoint,
	snvSiteSnapPrimerRefPoint,
	primerLensSnapPrimerRefPoint,
	snapbackBaseAtSNV,
	matchesWild,
	targetSnapMeltTemp
) {
	//──────────────────────────────────────────────────────────────────────────//
	//							Parameter Checking								//
	//──────────────────────────────────────────────────────────────────────────//
	// 1. targetStrandSeqSnapPrimerRefPoint
	if (!isValidDNASequence(targetStrandSeqSnapPrimerRefPoint)) {
		throw new Error(
			`Invalid targetStrandSeqSnapPrimerRefPoint: "${targetStrandSeqSnapPrimerRefPoint}". Must be a valid DNA sequence.`
		);
	}

	// 2. snvSiteSnapPrimerRefPoint
	if (
		typeof snvSiteSnapPrimerRefPoint !== 'object' ||
		snvSiteSnapPrimerRefPoint == null ||
		Array.isArray(snvSiteSnapPrimerRefPoint)
	) {
		throw new Error(
			`snvSiteSnapPrimerRefPoint must be a non-null object with "index" and "variantBase" keys.`
		);
	}
	const snvKeys = Object.keys(snvSiteSnapPrimerRefPoint);
	const validSNVKeys = new Set(['index', 'variantBase']);
	for (const key of snvKeys) {
		if (!validSNVKeys.has(key)) {
			throw new Error(
				`Unexpected key in snvSiteSnapPrimerRefPoint: "${key}"`
			);
		}
	}
	for (const expectedKey of validSNVKeys) {
		if (!(expectedKey in snvSiteSnapPrimerRefPoint)) {
			throw new Error(
				`Missing key in snvSiteSnapPrimerRefPoint: "${expectedKey}"`
			);
		}
	}
	if (
		typeof snvSiteSnapPrimerRefPoint.index !== 'number' ||
		!Number.isInteger(snvSiteSnapPrimerRefPoint.index) ||
		snvSiteSnapPrimerRefPoint.index < 0 ||
		snvSiteSnapPrimerRefPoint.index >=
			targetStrandSeqSnapPrimerRefPoint.length
	) {
		throw new Error(
			`snvSiteSnapPrimerRefPoint.index must be an integer in range [0, ${
				targetStrandSeqSnapPrimerRefPoint.length - 1
			}].`
		);
	}
	if (
		typeof snvSiteSnapPrimerRefPoint.variantBase !== 'string' ||
		snvSiteSnapPrimerRefPoint.variantBase.length !== 1 ||
		!VALID_BASES.has(snvSiteSnapPrimerRefPoint.variantBase)
	) {
		throw new Error(
			`snvSiteSnapPrimerRefPoint.variantBase must be one of "A", "T", "C", or "G".`
		);
	}

	// 3. snapbackBaseAtSNV
	if (
		typeof snapbackBaseAtSNV !== 'string' ||
		snapbackBaseAtSNV.length !== 1 ||
		!VALID_BASES.has(snapbackBaseAtSNV)
	) {
		throw new Error(
			`snapbackBaseAtSNV must be a single character base: "A", "T", "C", or "G".`
		);
	}

	// 4. matchesWild
	if (typeof matchesWild !== 'boolean') {
		throw new Error(`matchesWild must be a boolean.`);
	}

	// 5. primerLensSnapPrimerRefPoint
	if (
		typeof primerLensSnapPrimerRefPoint !== 'object' ||
		primerLensSnapPrimerRefPoint == null ||
		Array.isArray(primerLensSnapPrimerRefPoint)
	) {
		throw new Error(
			`primerLensSnapPrimerRefPoint must be an object with "primerLen" and "compPrimerLen" keys.`
		);
	}
	const primerKeys = Object.keys(primerLensSnapPrimerRefPoint);
	const validPrimerKeys = new Set(['primerLen', 'compPrimerLen']);
	for (const key of primerKeys) {
		if (!validPrimerKeys.has(key)) {
			throw new Error(
				`Unexpected key in primerLensSnapPrimerRefPoint: "${key}"`
			);
		}
	}
	for (const expectedKey of validPrimerKeys) {
		if (!(expectedKey in primerLensSnapPrimerRefPoint)) {
			throw new Error(
				`Missing key in primerLensSnapPrimerRefPoint: "${expectedKey}"`
			);
		}
	}
	for (const [key, val] of Object.entries(primerLensSnapPrimerRefPoint)) {
		if (
			typeof val !== 'number' ||
			!Number.isInteger(val) ||
			val < MIN_PRIMER_LEN
		) {
			throw new Error(
				`${key} must be an integer >= MIN_PRIMER_LEN (${MIN_PRIMER_LEN}). Got: ${val}`
			);
		}
	}

	// 6. SNV must be at least SNV_BASE_BUFFER away from both primers
	const { index: snvIndex } = snvSiteSnapPrimerRefPoint;
	const { primerLen, compPrimerLen } = primerLensSnapPrimerRefPoint;
	const seqLen = targetStrandSeqSnapPrimerRefPoint.length;

	if (
		snvIndex < primerLen + SNV_BASE_BUFFER ||
		snvIndex > seqLen - compPrimerLen - SNV_BASE_BUFFER - 1
	) {
		throw new Error(
			`SNV at index ${snvIndex} is too close to one of the primers. Must be at least ${SNV_BASE_BUFFER} bases away from both.`
		);
	}

	// 7. targetSnapMeltTemp
	if (
		typeof targetSnapMeltTemp !== 'number' ||
		!Number.isFinite(targetSnapMeltTemp) ||
		targetSnapMeltTemp <= 0
	) {
		throw new Error(
			`targetSnapMeltTemp must be a positive finite number. Got: ${targetSnapMeltTemp}`
		);
	}

	//──────────────────────────────────────────────────────────────────────────//
	//								Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//

	//// We will keep enlarging the stem, right then left... (as long as we are not up against the primers), keeping
	//// track of the melting temperature of the snapback for the wild-type allele, until we go over the desired meling
	//// temperature or we run out of viable stem location

	// Initialize variable to hold the snapback melting temperature for the wild and variant type alleles that is closest to the desired
	// snapback melting temperature for the wild type allele
	let bestWildTm = null;
	let bestVariantTm = null;
	// Initialize the variable for the corresponding stem locations
	let bestStemLocation = { start: null, end: null };

	// Initialize stem region
	let stemStart = snvIndex - SNV_BASE_BUFFER;
	let stemEnd = snvIndex + SNV_BASE_BUFFER;

	// Loop to grow stem until we go above desired Tm OR we've come up against both primers
	while (true) {
		// Slice current stem
		const currentStem = seq.slice(stemStart, stemEnd + 1);
		// Calculating the loop length
		const loopLen =
			stemStart + INNER_LOOP_NUMBER_OF_STRONG_BASE_MISMATCHES_REQUIRED;

		// Build mismatch object for wild and variant type, if needed, for stem Tm calculation
		let wildMismatch = null;
		let variantMismatch = null;
		if (!matchesWild) {
			wildMismatch = {
				position: snvIndex - stemStart, // Appropriate position is relative to the start of the stem
				type: snapbackBaseAtSNV,
			};
		} else {
			variantMismatch = {
				position: snvIndex - stemStart, // Appropriate position is relative to the start of the stem
				type: snapbackBaseAtSNV,
			};
		}

		// Compute the wild type allele melting temperature of the snapback
		const wildTm = await calculateSnapbackTm(
			currentStem,
			loopLen,
			wildMismatch
		);

		// Update the closest to desired wild type snapback melting temperature and corresponding stem location if applicable
		if (
			!bestWildTm ||
			Math.abs(wildTm - targetSnapMeltTemp) <
				Math.abs(bestWildTm - targetSnapMeltTemp)
		) {
			bestWildTm = wildTm;
			bestStemLocation.start = stemStart;
			bestStemLocation.end = stemEnd;
			// Compute and save the corresponding best variant type snapback temperature
			variantStemTm = await calculateSnapbackTm(
				currentStem,
				loopLen,
				variantMismatch
			);
		}

		// Loop Temination if wildTm has become larger than the desired snapback melting temperature
		// for the wild type allele
		if (wildTm >= targetSnapMeltTemp) {
			break;
		}

		// Grow the stem in the appropriate direction (if it can be grown without overlapping a primer location)
		if (
			snvIndex - stemStart < stemEnd - snvIndex &&
			stemStart > primerLen
		) {
			// We should push the start of the stem one nucleotide to the left
			stemStart -= 1;
		} else if (
			stemEnd - snvIndex <= snvIndex - stemStart &&
			stemEnd < seqLen - compPrimerLen - 1
		) {
			// We should push the start of the stem one nucleotide to the left
			stemEnd += 1;
		} else {
			// stem is up against both primers and we need to terminate the loop
			break;
		}
	}

	// Final check if final stem doesn’t meet minimum melting temperature requirement
	if (finalWildTm < MINIMUM_TARGET_SNAPBACK_MELTING_TEMP) {
		throw new Error(
			`Could not meet minimum snapback melting temp of ${MINIMUM_TARGET_SNAPBACK_MELTING_TEMP}°C. Final wildTm = ${finalWildTm.toFixed(
				2
			)}°C. Please consider moving primers farther out so a larger, more stable snapback stem can be created. `
		);
	}

	return {
		bestStemLocation,
		meltingTemps: {
			wildTm: parseFloat(bestWildTm.toFixed(TM_DECIMAL_PLACES)),
			variantTm: parseFloat(finalVariantTm.toFixed(TM_DECIMAL_PLACES)),
		},
		snapbackBaseAtSNV,
	};
}

/**
 * Produces the final snapback primer by:
 *   1. Taking the primer region on the chosen strand.
 *   2. Appending the reverse-complement of the stem (with the designated
 *      snapback base at the SNV position).
 *   3. Reversing the entire string so the resulting sequence is written 3'→5'.
 *
 * The returned string is exactly what is synthesized: snapback tail (3' end)
 * followed by the primer (5' end when annealed to the template).
 *
 * @typedef {Object} SNVSiteRefPoint
 * @property {number} index        0-based SNV index
 * @property {string} variantBase  Variant base at that position ("A", "T", "C", or "G").
 *
 * @typedef {Object} PrimerLensRefPoint
 * @property {number} primerLen        Length of the primer on this strand.
 * @property {number} compPrimerLen    Length of the complementary/limiting primer.
 *
 * @typedef {Object} StemLoc
 * @property {number} start  Inclusive start index of the stem.
 * @property {number} end    Inclusive end index of the stem.
 *
 * @param {string} targetStrandSeqSnapPrimerRefPoint - DNA sequence (5'→3') for the strand that receives the snapback.
 * @param {SNVSiteRefPoint} snvSiteSnapPrimerRefPoint - SNV description in reference to the target strand
 * @param {PrimerLensRefPoint} primerLensSnapPrimerRefPoint - Object holding the primer length and the limiting primer length
 * @param {StemLoc} stemLoc - Location of the stem on this strand.
 * @param {string} snapbackBaseAtSNV - Base placed in the snapback tail at the SNV site.
 *
 * @returns {string}  Final snapback primer sequence written 3'→5'.
 */
function buildFinalSnapback(
	targetStrandSeqSnapPrimerRefPoint,
	snvSiteSnapPrimerRefPoint,
	primerLensSnapPrimerRefPoint,
	stemLoc,
	snapbackBaseAtSNV
) {
	//──────────────────────────────────────────────────────────────────────────//
	//							Parameter Checking								//
	//──────────────────────────────────────────────────────────────────────────//
	// 1. targetStrandSeqSnapPrimerRefPoint
	if (!isValidDNASequence(targetStrandSeqSnapPrimerRefPoint)) {
		throw new Error(
			`Invalid targetStrandSeqSnapPrimerRefPoint: "${targetStrandSeqSnapPrimerRefPoint}". Must be a valid DNA sequence.`
		);
	}

	// 2. snvSiteSnapPrimerRefPoint
	if (
		typeof snvSiteSnapPrimerRefPoint !== 'object' ||
		snvSiteSnapPrimerRefPoint == null ||
		Array.isArray(snvSiteSnapPrimerRefPoint)
	) {
		throw new Error(
			`snvSiteSnapPrimerRefPoint must be a non-null object with "index" and "variantBase" keys.`
		);
	}
	const snvKeys = Object.keys(snvSiteSnapPrimerRefPoint);
	const validSNVKeys = new Set(['index', 'variantBase']);
	for (const key of snvKeys) {
		if (!validSNVKeys.has(key)) {
			throw new Error(
				`Unexpected key in snvSiteSnapPrimerRefPoint: "${key}"`
			);
		}
	}
	for (const expectedKey of validSNVKeys) {
		if (!(expectedKey in snvSiteSnapPrimerRefPoint)) {
			throw new Error(
				`Missing key in snvSiteSnapPrimerRefPoint: "${expectedKey}"`
			);
		}
	}
	if (
		typeof snvSiteSnapPrimerRefPoint.index !== 'number' ||
		!Number.isInteger(snvSiteSnapPrimerRefPoint.index) ||
		snvSiteSnapPrimerRefPoint.index < 0 ||
		snvSiteSnapPrimerRefPoint.index >=
			targetStrandSeqSnapPrimerRefPoint.length
	) {
		throw new Error(
			`snvSiteSnapPrimerRefPoint.index must be an integer in range [0, ${
				targetStrandSeqSnapPrimerRefPoint.length - 1
			}].`
		);
	}
	if (
		typeof snvSiteSnapPrimerRefPoint.variantBase !== 'string' ||
		snvSiteSnapPrimerRefPoint.variantBase.length !== 1 ||
		!VALID_BASES.has(snvSiteSnapPrimerRefPoint.variantBase)
	) {
		throw new Error(
			`snvSiteSnapPrimerRefPoint.variantBase must be one of "A", "T", "C", or "G".`
		);
	}

	// 3. snapbackBaseAtSNV
	if (
		typeof snapbackBaseAtSNV !== 'string' ||
		snapbackBaseAtSNV.length !== 1 ||
		!VALID_BASES.has(snapbackBaseAtSNV)
	) {
		throw new Error(
			`snapbackBaseAtSNV must be a single character base: "A", "T", "C", or "G".`
		);
	}

	// 5. primerLensSnapPrimerRefPoint
	if (
		typeof primerLensSnapPrimerRefPoint !== 'object' ||
		primerLensSnapPrimerRefPoint == null ||
		Array.isArray(primerLensSnapPrimerRefPoint)
	) {
		throw new Error(
			`primerLensSnapPrimerRefPoint must be an object with "primerLen" and "compPrimerLen" keys.`
		);
	}
	const primerKeys = Object.keys(primerLensSnapPrimerRefPoint);
	const validPrimerKeys = new Set(['primerLen', 'compPrimerLen']);
	for (const key of primerKeys) {
		if (!validPrimerKeys.has(key)) {
			throw new Error(
				`Unexpected key in primerLensSnapPrimerRefPoint: "${key}"`
			);
		}
	}
	for (const expectedKey of validPrimerKeys) {
		if (!(expectedKey in primerLensSnapPrimerRefPoint)) {
			throw new Error(
				`Missing key in primerLensSnapPrimerRefPoint: "${expectedKey}"`
			);
		}
	}
	for (const [key, val] of Object.entries(primerLensSnapPrimerRefPoint)) {
		if (
			typeof val !== 'number' ||
			!Number.isInteger(val) ||
			val < MIN_PRIMER_LEN
		) {
			throw new Error(
				`${key} must be an integer >= MIN_PRIMER_LEN (${MIN_PRIMER_LEN}). Got: ${val}`
			);
		}
	}

	// 6) stemLoc
	if (
		typeof stemLoc !== 'object' ||
		stemLoc == null ||
		Array.isArray(stemLoc) ||
		!('start' in stemLoc) ||
		!('end' in stemLoc)
	) {
		throw new Error(
			`stemLoc must be an object with "start" and "end" properties.`
		);
	}
	{
		const validStemKeys = new Set(['start', 'end']);
		for (const key of Object.keys(stemLoc)) {
			if (!validStemKeys.has(key)) {
				throw new Error(`Unexpected key in stemLoc: "${key}"`);
			}
		}
	}
	if (
		typeof stemLoc.start !== 'number' ||
		typeof stemLoc.end !== 'number' ||
		!Number.isInteger(stemLoc.start) ||
		!Number.isInteger(stemLoc.end) ||
		stemLoc.start < 0 ||
		stemLoc.end < 0 ||
		stemLoc.start >= targetStrandSeqSnapPrimerRefPoint.length ||
		stemLoc.end >= targetStrandSeqSnapPrimerRefPoint.length ||
		stemLoc.start > stemLoc.end
	) {
		throw new Error(
			`stemLoc.start and stemLoc.end must be valid indices within sequence bounds and start ≤ end.`
		);
	}

	// 7) make sure we can still add the required strong-mismatch bases
	//    without falling off the end of the sequence
	const maxIndexNeeded =
		stemLoc.end + END_OF_STEM_NUMBER_OF_STRONG_BASE_MISMATCHES_REQUIRED;

	if (maxIndexNeeded >= targetStrandSeqSnapPrimerRefPoint.length) {
		throw new Error(
			`stemLoc.end (${stemLoc.end}) + ${END_OF_STEM_NUMBER_OF_STRONG_BASE_MISMATCHES_REQUIRED} strong-mismatch bases exceeds sequence length ` +
				`(${targetStrandSeqSnapPrimerRefPoint.length}).`
		);
	}

	//──────────────────────────────────────────────────────────────────────────//
	//								Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//
	// aliases
	const seq = targetStrandSeqSnapPrimerRefPoint;
	const primerLen = primerLensSnapPrimerRefPoint.primerLen;
	const snvIndex = snvSiteSnapPrimerRefPoint.index;
	const { stemStart, stemEnd } = stemLoc;

	// 1) Start with the primer and add tail to the left of it
	let snapback = seq.slice(0, primerLen);

	// 2) Add the inner loop strong mismatches to the left of the snapback
	for (
		let i =
			stemStart - INNER_LOOP_NUMBER_OF_STRONG_BASE_MISMATCHES_REQUIRED;
		i < stemStart;
		i--
	) {
		snapback = STRONG_NUCLEOTIDE_MISMATCH[seq[i]] + snapback;
	}

	// 3) Add the stem sequence, from the start to the end, to the beggining of the snapback.
	//    If the site is the SNV site we add the snapback base chosen for the SNV site
	for (let i = stemStart; i <= stemEnd; i++) {
		if (i === snvIndex) {
			snapback = snapbackBaseAtSNV + snapback;
		} else {
			snapback = NUCLEOTIDE_COMPLEMENT[seq[i]] + snapback;
		}
	}

	// 4) Add the end of stem strong mismatched
	for (
		let i = stemEnd + 1;
		i <= stemEnd + END_OF_STEM_NUMBER_OF_STRONG_BASE_MISMATCHES_REQUIRED;
		i++
	) {
		snapback = STRONG_NUCLEOTIDE_MISMATCH[seq[i]] + snapback;
	}

	return snapback;
}

/*****************************************************************************************/
/************************************ Helper Function ************************************/
/*****************************************************************************************/

/**
 * Checks if the SNV is too close to the ends of the primers to form 3 base buffer on either end.
 *
 * Assumptions:
 * - The target sequence strand starts at the 5' end
 * - Primer corresponds to primer that binds to the target sequence strand given
 * - Complementary primer corresponds to the complement of target sequence strand given
 *
 * @param {number} snvIndex - Index of the SNV in the target sequence strand(starting at 0).
 * @param {number} primerLen - Length of the primer on the same strand as the target sequence.
 * @param {number} compPrimerLen - Length of the complementary primer on the opposite strand.
 * @param {number} seqLen - The total length of the target sequence.
 *
 * @returns {boolean} - True if SNV is within 3 bases of either primer (i.e. too close), otherwise false.
 * @throws {Error} - If any argument is missing, invalid, or out of bounds.
 */
function snvTooCloseToPrimer(snvIndex, primerLen, compPrimerLen, seqLen) {
	//──────────────────────────────────────────────────────────────────────────//
	//							Parameter Checking								//
	//──────────────────────────────────────────────────────────────────────────//
	// Check for missing or non-numeric inputs
	for (const [name, val] of [
		['snvIndex', snvIndex],
		['primerLen', primerLen],
		['compPrimerLen', compPrimerLen],
		['seqLen', seqLen],
	]) {
		if (typeof val !== 'number' || !Number.isFinite(val)) {
			throw new Error(`${name} must be a finite number`);
		}
		if (!Number.isInteger(val)) {
			throw new Error(`${name} must be an integer`);
		}
		if (val < 0) {
			throw new Error(`${name} must be non-negative`);
		}
	}
	// Validate SNV index is in bounds
	if (snvIndex >= seqLen) {
		throw new Error(
			`snvIndex (${snvIndex}) cannot exceed or equal sequence length (${
				seqLen - 1
			})`
		);
	}
	// Validate primer lengths are feasible (or there was a bigger mistake)
	if (primerLen + compPrimerLen >= seqLen) {
		throw new Error(
			`Primer lengths (${primerLen} + ${compPrimerLen}) exceed or match sequence length (${seqLen})`
		);
	}

	//──────────────────────────────────────────────────────────────────────────//
	//								Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//

	// Must be at least {SNV_BASE_BUFFER} bases away from either primer
	const lowerBoundIndex = primerLen + SNV_BASE_BUFFER;
	const upperBoundIndex = seqLen - compPrimerLen - SNV_BASE_BUFFER - 1;

	return snvIndex < lowerBoundIndex || snvIndex > upperBoundIndex;
}

/**
 * Constructs the mmseq string needed by the Tm service so it sees the intended
 * mismatch in the final double-stranded structure.
 *
 * For example, if mismatch.type = 'G', that means you want an A↔G mismatch in
 * the final pairing. The Tm service expects to see the difference as:
 *   seq=... 'A' ...
 *   mmseq=... 'C' ... (the complement of 'G') at that same position.
 *
 * @typedef {Object} Mismatch
 * @property {number} position - The index in `seq` where the mismatch occurs.
 * @property {string} type     - The base you'll have on the opposite strand (e.g. "G") (mismatch).
 *
 *
 * @param {string} seq - Original matched sequence (5'→3').
 * @param {Mismatch} mismatch - The mismatch specification.
 *
 * @throws {Error} If the sequence or mismatch is invalid.
 * @returns {string} The mmseq string for the Tm service.
 */
function buildMismatchSequenceForAPI(seq, mismatch) {
	//──────────────────────────────────────────────────────────────────────────//
	//							Parameter Checking								//
	//──────────────────────────────────────────────────────────────────────────//
	// Sequence must be a valid DNA string
	if (!isValidDNASequence(seq)) {
		throw new Error(
			`Invalid input sequence: "${seq}". Must only contain A, T, C, or G and be a string.`
		);
	}
	// Mismatch must be a javascript object with position and type keys
	if (
		typeof mismatch !== 'object' ||
		mismatch == null ||
		Array.isArray(mismatch) ||
		!('position' in mismatch) ||
		!('type' in mismatch)
	) {
		throw new Error(
			`Invalid mismatch object. Expected { position: number, type: string }. Received: ${JSON.stringify(
				mismatch
			)}`
		);
	}
	// mismatch.position must be an integer in [0, seqLen-1]
	if (
		typeof mismatch.position !== 'number' ||
		!Number.isInteger(mismatch.position) ||
		mismatch.position < 0 ||
		mismatch.position >= seq.length
	) {
		throw new Error(
			`Mismatch position ${mismatch.position} is not an integer or out of bounds for sequence of length ${seq.length}`
		);
	}
	// mismatch.type must be a string with one character: 'A', 'C', 'T', or 'G'
	if (
		typeof mismatch.type !== 'string' ||
		mismatch.type.length !== 1 ||
		!VALID_BASES.has(mismatch.type)
	) {
		throw new Error(
			`Mismatch type "${mismatch.type}" must be a single character: A, T, C, or G`
		);
	}

	//──────────────────────────────────────────────────────────────────────────//
	//								Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//

	// Get the complement base
	// Should never throw this error as parameters were already checked but can leave in, in case
	// code structure changes in the future
	const complementBase = NUCLEOTIDE_COMPLEMENT[mismatch.type];
	if (!complementBase) {
		throw new Error(
			`Unexpected error: unable to determine complement for mismatch base "${mismatch.type}". Error shouldn't be 
			thrown as parameters are checked.`
		);
	}

	// Replace that position with the mismatch complement's in the original sequence, this is what the API wants
	return (
		seq.slice(0, mismatch.position) +
		complementBase +
		seq.slice(mismatch.position + 1)
	);
}

/**
 * Simple parser to extract the numeric Tm from the raw HTML:
 * e.g. <html><head></head><body><seq>...</seq><tm>47.27</tm><mmtm>37.54</mmtm></body></html>
 *
 * @param {string} rawHtml - The raw string returned by .cgi file
 * @param {boolean} [mismatch] - Optional specification denoting if we want to parse out the mismatched tm
 * @returns {number|null}  - The Tm if found, otherwise null
 */
function parseTmFromResponse(rawHtml, mismatch) {
	//──────────────────────────────────────────────────────────────────────────//
	//							Parameter Checking								//
	//──────────────────────────────────────────────────────────────────────────//

	//──────────────────────────────────────────────────────────────────────────//
	//								Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//
	try {
		// Parse the text into a DOM (there might be a better way to do this with functionality built into Node too)
		const parser = new DOMParser();
		const doc = parser.parseFromString(rawHtml, 'text/html');

		// Getting the <tm> or <mmtm> element
		var tmElement;
		if (!mismatch) {
			// Look for a <tm> element
			tmElement = doc.querySelector('tm');
		} else {
			tmElement = doc.querySelector('mmtm');
		}

		// Returns null if correct tm is not found
		if (!tmElement) {
			return null;
		}

		// Convert the text inside element to a float
		const tmValue = parseFloat(tmElement.textContent.trim());
		return isNaN(tmValue) ? null : tmValue;
	} catch (err) {
		console.error('parseTmFromResponse error:', err);
		return null;
	}
}

/**
 * Estimates the melting temperature of a snapback structure using the formula:
 *
 *     Tm = -5.25 * ln(loopLen) + 0.837 * stemTm + 32.9
 *
 * stemTm is calculated via `getStemTm()`.
 *
 * Assumptions:
 * - The stem sequence is 5'→3' and corresponds to the strand containing the SNV, not the
 *   snapback tail strand.
 * - The mismatch position is within bounds of the stem.
 * - The snapback base is the one used in the snapback at the SNV position.
 * - The loopLen is at least {MIN_LOOP_LEN} bases long
 *
 * @typedef {Object} Mismatch
 * @property {number} position - Index in the stem where mismatch occurs.
 * @property {string} type     - Base on the opposite strand at the mismatch site
 * 								 This would be the snapback base at the SNV location, on the tail,
 * 								 if the wild type does not match the snapback tail.
 *
 * @param {string} stemSeq - The stem sequence (5'→3'), which includes the SNV.
 * @param {number} loopLen - The length (in nucleotides) of the snapback loop (unpaired region).
 * @param {Mismatch} [mismatch] - Optional mismatch object specifying the mismatch position and base.
 *
 * @returns {Promise<number>} - The estimated melting temperature (Tm) in degrees Celsius.
 *
 * @throws {Error} - If any input is missing or invalid.
 */
async function calculateSnapbackTm(stemSeq, loopLen, mismatch) {
	//──────────────────────────────────────────────────────────────────────────//
	//							Parameter Checking								//
	//──────────────────────────────────────────────────────────────────────────//
	// Checking the stemSeq parameter
	if (!isValidDNASequence(stemSeq)) {
		throw new Error(`Invalid DNA sequence passed to stemSeq: "${stemSeq}"`);
	}
	// Checking the loopLen parameter
	if (
		typeof loopLen !== 'number' ||
		!Number.isFinite(loopLen) ||
		loopLen < MIN_LOOP_LEN
	) {
		throw new Error(
			`loopLen must be a positive, finite number greater than or equal to ${MIN_LOOP_LEN}`
		);
	}
	// Checking the mismatch parameter, if it is passed
	if (mismatch !== undefined && mismatch !== null) {
		// Check its a JavaScript object
		if (typeof mismatch !== 'object' || Array.isArray(mismatch)) {
			throw new Error(
				`Mismatch must be an object if provided. Received: ${typeof mismatch}`
			);
		}
		// Check that it contains the right keys
		const allowedKeys = new Set(['position', 'type']);
		const actualKeys = Object.keys(mismatch);
		for (const key of actualKeys) {
			if (!allowedKeys.has(key)) {
				throw new Error(
					`Mismatch object contains unexpected key: "${key}"`
				);
			}
		}
		if (!('position' in mismatch) || !('type' in mismatch)) {
			throw new Error(
				`Mismatch object must include both "position" and "type". Got: ${JSON.stringify(
					mismatch
				)}`
			);
		}
		//// Check the key values
		// Validate position type
		if (
			typeof mismatch.position !== 'number' ||
			!Number.isInteger(mismatch.position)
		) {
			throw new Error(`Mismatch.position must be an integer.`);
		}
		// Ensure mismatch is not too close to ends of the stem
		if (
			mismatch.position < SNV_BASE_BUFFER ||
			mismatch.position > stemSeq.length - SNV_BASE_BUFFER - 1
		) {
			throw new Error(
				`Mismatch.position (${mismatch.position}) must be at least ${SNV_BASE_BUFFER} bases from both ends of the stem (length: ${stemSeq.length})`
			);
		}
		// Validate mismatch type is a nucleotide
		if (
			typeof mismatch.type !== 'string' ||
			mismatch.type.length !== 1 ||
			!VALID_BASES.has(mismatch.type)
		) {
			throw new Error(
				`Mismatch.type must be a single character: A, T, C, or G. Got: "${mismatch.type}"`
			);
		}
	}

	//──────────────────────────────────────────────────────────────────────────//
	//								Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//

	// Calculate the wild type's stem melting temperature, only pass the mismatch if it was passed
	const stemTM = await getStemTm(stemSeq, mismatch ?? undefined);

	// Plug into the Tm formula
	const tm = -5.25 * Math.log(loopLen) + 0.837 * stemTM + 32.9;

	// Return rounded to {TM_DECIMAL_PLACES} decimal place(s)
	return parseFloat(tm.toFixed(TM_DECIMAL_PLACES));
}

/*****************************************************************************************/
/********************************** DNA Utility Function *********************************/
/*****************************************************************************************/

/**
 * Determines whether `seqStrand` is a valid, non-empty, uppercase
 * DNA string consisting solely of the characters A, T, C, or G.
 *
 * ─────────────────────────────────────────────────────────────────────────────
 * Assumptions
 * ─────────────────────────────────────────────────────────────────────────────
 * - A valid DNA sequence must be uppercase and contain at least one nucleotide.
 *
 * ─────────────────────────────────────────────────────────────────────────────
 * Parameters, Returns, and Errors
 * ─────────────────────────────────────────────────────────────────────────────
 * @param   {string}  seqStrand   Candidate DNA sequence to validate.
 *
 * @returns {boolean}             true  → valid DNA string
 *                                false → invalid type, empty, or contains
 *                                        characters outside A/T/C/G.
 */
function isValidDNASequence(seqStrand) {
	//──────────────────────────────────────────────────────────────────────//
	//                           Function Logic                             //
	//──────────────────────────────────────────────────────────────────────//
	if (typeof seqStrand !== 'string' || seqStrand.length === 0) {
		return false;
	}

	for (const base of seqStrand) {
		if (!VALID_BASES.has(base)) {
			return false;
		}
	}

	return true;
}

/**
 * Returns the complement of a DNA sequence.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Assumptions
 * ──────────────────────────────────────────────────────────────────────────
 * - Input must be a valid uppercase DNA sequence consisting of characters:
 *   A, T, C, or G.
 * - The complement sequence is NOT reversed. For example, the complement of 'GA'
 *   is 'CT'.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Parameters, Returns, and Errors
 * ──────────────────────────────────────────────────────────────────────────
 * @param   {string}   seqStrand     DNA sequence (e.g., "ATCG").
 *
 * @returns {string}                 Complement of the input (e.g., "TAGC").
 *
 * @throws  {Error}                  If input is not a valid DNA sequence.
 */
function complementSequence(seqStrand) {
	//──────────────────────────────────────────────────────────────────────────//
	//							Parameter Checking								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1. seqStrand
	if (!isValidDNASequence(seqStrand)) {
		throw new Error(
			`Invalid DNA sequence: ${seqStrand}. ` +
				`Must be a non-empty uppercase string containing only characters A, T, C, and/or G.`
		);
	}

	//──────────────────────────────────────────────────────────────────────────//
	//								Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1. Convert each base to its complement
	return seqStrand
		.split('')
		.map((base) => NUCLEOTIDE_COMPLEMENT[base])
		.join('');
}

/**
 * Returns the reverse complement of a DNA sequence.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Assumptions
 * ──────────────────────────────────────────────────────────────────────────
 * - The input sequence is a valid, uppercase DNA sequence using only A, T, C, or G.
 * - The reverse complement is defined as:
 *     1. Replacing each base with its complement (A<->T, C<->G)
 *     2. Reversing the entire resulting string to keep 5' to 3' orientation
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Parameters, Returns, and Errors
 * ──────────────────────────────────────────────────────────────────────────
 * @param   {string}   seqStrand     DNA sequence (e.g., "ATCG").
 *
 * @returns {string}                 Reverse complement of `seqStrand` (e.g., "CGAT").
 *
 * @throws  {Error}                  If input is empty, not a string, or contains
 *                                   invalid DNA characters.
 */
function reverseComplement(seqStrand) {
	//──────────────────────────────────────────────────────────────────────────//
	//							Parameter Checking								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1. seqStrand
	if (!isValidDNASequence(seqStrand)) {
		throw new Error(
			`Invalid DNA sequence: ${seqStrand}. ` +
				`Must be a non-empty uppercase string containing only characters A, T, C, and/or G.`
		);
	}

	//──────────────────────────────────────────────────────────────────────────//
	//								Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1. Get complement of sequence
	const complementStrand = complementSequence(seqStrand);

	// 2. Reverse the complement strand
	return complementStrand.split('').reverse().join('');
}

/**
 * Determines whether a given object is a valid SNVSite object.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Assumptions
 * ──────────────────────────────────────────────────────────────────────────
 * - An SNV object must have exactly two keys:
 *     1. `index`: a non-negative integer (0 or greater)
 *     2. `variantBase`: one of the characters "A", "T", "C", or "G"
 * - The object must not contain any other keys.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Parameters, Returns, and Errors
 * ──────────────────────────────────────────────────────────────────────────
 * @param   {Object}  snv    Object to validate as an SNVSite.
 *
 * @returns {boolean}        true  → if object is a valid SNVSite
 *                           false → otherwise
 */
function isValidSNVObject(snv) {
	//──────────────────────────────────────────────────────────────────────────//
	//							Parameter Checking								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1. Type check: Must be a non-null object (and not an array)
	if (typeof snv !== 'object' || snv === null || Array.isArray(snv)) {
		return false;
	}

	// 2. Key check: Must contain exactly 'index' and 'variantBase'
	const expectedKeys = new Set(['index', 'variantBase']);
	const actualKeys = Object.keys(snv);

	if (actualKeys.length !== expectedKeys.size) {
		return false;
	}
	for (const key of actualKeys) {
		if (!expectedKeys.has(key)) {
			return false;
		}
	}

	// 3. Validate `index` is a non-negative integer
	if (
		typeof snv.index !== 'number' ||
		!Number.isInteger(snv.index) ||
		snv.index < 0
	) {
		return false;
	}

	// 4. Validate `variantBase` is a valid uppercase base
	if (
		typeof snv.variantBase !== 'string' ||
		snv.variantBase.length !== 1 ||
		!VALID_BASES.has(snv.variantBase)
	) {
		return false;
	}

	//──────────────────────────────────────────────────────────────────────────//
	//								Validation Passed							//
	//──────────────────────────────────────────────────────────────────────────//

	return true;
}

/**
 * Returns the reverse complement SNV site.
 * That is, it transforms the SNV's index and base to be correct
 * for the reverse complement strand.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Assumptions
 * ──────────────────────────────────────────────────────────────────────────
 * - `snvSite` must be a valid SNVSite object, verified by isValidSNVObject().
 * - `seqLen` must be a positive integer ≥ 1.
 * - `snvSite.index` must be in the range [0, seqLen - 1].
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Parameters, Returns, and Errors
 * ──────────────────────────────────────────────────────────────────────────
 * @param   {SNVSite} snvSite     Object with properties:
 *                                { index: number, variantBase: "A"|"T"|"C"|"G" }
 *
 * @param   {number}  seqLen      Length of the DNA sequence the SNV belongs to.
 *
 * @returns {SNVSite}             SNV object mapped to the reverse complement strand:
 *                                {
 *                                  index: number,
 *                                  variantBase: string
 *                                }
 *
 * @throws  {Error}               If snvSite is invalid or index is out of bounds,
 *                                or if seqLen is not a valid positive integer.
 */
function revCompSNV(snvSite, seqLen) {
	//──────────────────────────────────────────────────────────────────────────//
	//							Parameter Checking								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1. Validate SNV object structure
	if (!isValidSNVObject(snvSite)) {
		throw new Error(
			`Invalid SNV object: must be { index: number, variantBase: "A"|"T"|"C"|"G" }. Received: ${JSON.stringify(
				snvSite
			)}`
		);
	}

	// 2. Validate seqLen
	if (
		typeof seqLen !== 'number' ||
		!Number.isInteger(seqLen) ||
		seqLen <= 0
	) {
		throw new Error(
			`seqLen must be a positive integer. Received: ${seqLen}`
		);
	}

	// 3. Ensure index is within sequence bounds
	if (snvSite.index >= seqLen) {
		throw new Error(
			`snvSite.index (${snvSite.index}) is out of bounds for sequence length ${seqLen}`
		);
	}

	//──────────────────────────────────────────────────────────────────────────//
	//								Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1. Compute reverse complement index
	const revCompIndex = seqLen - snvSite.index - 1;

	// 2. Convert variant base to its complement
	const revCompBase = NUCLEOTIDE_COMPLEMENT[snvSite.variantBase];

	// 3. Return updated SNV object
	return {
		index: revCompIndex,
		variantBase: revCompBase,
	};
}

/**
 * Reverses a DNA sequence
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Assumptions
 * ──────────────────────────────────────────────────────────────────────────
 * - `seqStrand` is expected to be a non-empty uppercase DNA string containing
 *   only the characters “A”, “T”, “C”, or “G”.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Parameters, Returns, and Errors
 * ──────────────────────────────────────────────────────────────────────────
 * @param   {string}  		seqStrand   	DNA sequence to reverse.
 *
 * @returns {string}              			The reversed sequence.
 *
 * @throws  {Error}               			If `seqStrand` is empty, not a string,
 *                                			or contains characters other than A/T/C/G.
 */
function reverseSequence(seqStrand) {
	//──────────────────────────────────────────────────────────────────────────//
	//							Parameter Checking								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1. seqStrand
	if (!isValidDNASequence(seqStrand)) {
		throw new Error(
			`Invalid DNA sequence: ${seqStrand}. ` +
				`Must be a non-empty uppercase string containing only characters A, T, C, and/or G..`
		);
	}

	//──────────────────────────────────────────────────────────────────────────//
	//								Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1. Reverse the sequence
	return seqStrand.split('').reverse().join('');
}

/*****************************************************************************************/
/************************************ Export Function ************************************/
/*****************************************************************************************/

// For testing with script.test.js
export {
	// Primary function
	createSnapback,

	// Secondary functions
	useTargetStrandsPrimerForComplement,
	evaluateSnapbackMatchingOptions,
	getStemTm,
	createStem,
	buildFinalSnapback,

	// Helper/logic functions
	snvTooCloseToPrimer,
	buildMismatchSequenceForAPI,
	parseTmFromResponse,
	calculateSnapbackTm,

	// DNA utility functions
	isValidDNASequence,
	isValidSNVObject,
	complementSequence,
	reverseComplement,
	revCompSNV,
	reverseSequence,

	// Constants
	SNV_BASE_BUFFER,
	NUCLEOTIDE_COMPLEMENT,
	MIN_LOOP_LEN,
	MIN_PRIMER_LEN,
	END_OF_STEM_NUMBER_OF_STRONG_BASE_MISMATCHES_REQUIRED,
};
