/*
File:           script.js
Description:    Main JavaScript file for interactivity, functionality, and processing.
Author:         Adrian deCola
Relative Path:  uSnapback/src/script.js
*/

/****************************************************************/
/************************** Constants ***************************/
/****************************************************************/
const NUCLEOTIDE_COMPLEMENT = { A: 'T', T: 'A', C: 'G', G: 'C' };
const STRONG_NUCLEOTIDE_MISMATCH = { A: 'G', G: 'A', C: 'C', T: 'T' };
const VALID_BASES = new Set(['A', 'T', 'C', 'G']);
const SNV_BASE_BUFFER = 4; // The number of matched bases required on either end of a mismatched SNV
const INNER_LOOP_NUMBER_OF_STRONG_BASE_MISMATCHES_REQUIRED = 2;
const END_OF_STEM_NUMBER_OF_STRONG_BASE_MISMATCHES_REQUIRED = 2;
const MINIMUM_TARGET_SNAPBACK_MELTING_TEMP = 40;
const MAXIMUM_TARGET_SNAPBACK_MELTING_TEMP = 80;
const SNAP_DECIMAL_PLACES = 2;
const API_TOKEN = 1;
// Chemisty parameters ************** name these better ************
const MG = 2.2;
const MONO = 20.0;
const T_PARAM = 'UnifiedSantaLucia';
const SALT_CALC_TYPE = 'bpdenominator';
const O_TYPE = 'oligo';
const CONC = 0.5;
const LIMITING_CONC = 0.5;

/****************************************************************/
/*********************** Primary Function ***********************/
/****************************************************************/
/**
 * Creates a snapback primer sequence by identifying suitable stem regions
 * in the target DNA sequence, accounting for primer lengths and a single
 * mismatch site. The function then chooses which primer to add the snapback
 * tail to and whether it should match with the wild type of one of the
 * variants, if there are multiple. It then creates the rest of the stem
 * so that it gets as close as it can to the desired snapback melting temperature
 * and returns an error if it cannot have the snapback melting temperatures
 * (match or mismatch) all be aboge the minimum snapback melting temperature.
 *
 * Assumptions:
 * - Assumes passed in values are correct
 *     - Passed in target sequence is valid AND uppercase
 *
 * - Assumes that the target sequence begins with the 5' end
 * - we want to minimize loop length for whatever reason??
 * - primers should be at least 6 bases each (or at least the snapback one?)
 * - The stem's melting temperatuere is calcualted using the Santa Lucia method.
 * - Extension can occur on the hairpin's complement, on its 5' end. A 2 base pair mismatch after the stem can
 *   help avoid this
 *
 * Things to consider/add/make better/more effecient:
 *  * - I could represent DNA sequences as an array of 2 bit encoded neucleotides. This could save some memory
 * - I could save the melting temperatur so far and implement santa lucia method here and simply add entropys
 *   and enthalpies of added neucleotides.
 * - make my own data type to represent dna sequences with their own methods/properties like .complement or .complement()
 *
 *
 * @param {string} targetSeqStrand - A strand of the full DNA sequence to design the snapback primer for.
 * @param {number} primerLen - The length of the primer on the strand denoted by the target sequence
 * @param {number} compPrimerLen - The length of the complementary primer on the opposite strand.
 * @param {Object} snvSite - An object representing the single neucleotide variant site with:
 *     @property {number} index - The index in the sequence where the variant occurs (0-based).
 *     @property {string} variantBase - A one character string representing the variant base.
 * @param {number} minSnapbackMeltTemp - The minimum viable snapback melting temperature for any match or mismatch snapback
 * @param {number} desiredSnapbackMeltTempWildType = The desired snapback melting temperature for the wild type
 *
 * @returns {Object} - An object representing the formed snapback primer
 *     @property {string} snapbackSeq - The final snapback primer sequence (5'-tail + primer-3').
 *     @property {boolean} isOnTargetPrimer - A boolean representing if the primer to add a tail on is the primer denoted
 *                                            by the target sequence (True) or by its complement (False)
 *     @property {Array<Object>} snapbackMeltingTemps - Melting temperature info for match and mismatch variants.
 *         @property {string} snapbackMeltingTemps[].isWild - True if this is the melting temperature corresponding to the wild allele
 *         @property {boolean} snapbackMeltingTemps[].isMatch - True if this base is the correct one (no mismatch in snapback).
 *         @property {number} snapbackMeltingTemps[].meltingTemp - The melting temperature for this base configuration.
 *
 */
async function createSnapback(
	targetSeqStrand,
	primerLen,
	compPrimerLen,
	snvSite,
	desiredSnapbackMeltTempWildType
) {
	/////////// Parameter Checking ///////////
	// Validate targetSeqStrand is a valid DNA sequence
	if (!isValidDNASequence(targetSeqStrand)) {
		throw new Error(`Invalid targetSeqStrand: "${targetSeqStrand}"`);
	}

	// Validate desiredSnapbackMeltTempWildType is a positive number
	if (
		typeof desiredSnapbackMeltTempWildType !== 'number' ||
		!Number.isFinite(desiredSnapbackMeltTempWildType) ||
		desiredSnapbackMeltTempWildType < 0
	) {
		throw new Error(
			`desiredSnapbackMeltTempWildType must be a positive, finite number`
		);
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
	/////////// Funciton Logic ///////////

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
	const { useTargetStrand, snapbackBaseAtSNV } =
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
	const objQuestionMark = createStem(
		targetStrandSeqSnapPrimerRefPoint,
		snvSiteSnapPrimerRefPoint,
		snapbackBaseAtSNV,
		desiredSnapbackMeltTempWildType
	);

	// Pass in snapbackbase, stem locations, andINNER_LOOP_NUMBER_OF_STRONG_BASE_MISMATCHES_REQUIRED and END_OF_STEM_NUMBER_OF_STRONG_BASE_MISMATCHES_REQUIRED
	// to some build snapback function. get the snapback melting temperature, assumes nothing added to the end of stem or have it passed but createStem()
	//
}

/****************************************************************/
/********************* Secondary Functions **********************/
/****************************************************************/

/**
 * Decides which strand (target vs. complementary) should get the snapback primer,
 * and whether the snapback should match the wild-type base or the single variant base at the SNV.
 *
 * Returns an object:
 *   {
 *     useTargetStrand: boolean,   // True => use target strand for snapback
 *     snapbackBaseAtSNV: string,  // e.g. 'A' or 'G', whichever base is best to match
 *   }
 *
 * @param {string} targetSeqStrand     - Full target strand (5'→3')
 * @param {Object} snvSite             - { index: number, variantBase: string }
 *
 * @returns {Promise<{ useTargetStrand: boolean, snapbackBaseAtSNV: string }>}
 */
async function useTargetStrandsPrimerForComplement(targetSeqStrand, snvSite) {
	/////////// Parameter Checking ///////////

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

	/////////// Function Logic ///////////

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
		};
	} else {
		return {
			useTargetStrand: false,
			snapbackBaseAtSNV: compScenario.bestSnapbackBase,
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
 * @returns {Promise<{ bestSnapbackBase: string, bestDifference: number }>}
 *   bestSnapbackBase => the base to use in the snapback (either complementary to the wild or variant base)
 *   bestDifference => the numeric Tm difference from the wild scenario
 */
async function evaluateSnapbackMatchingOptions(
	initStem,
	mismatchPos,
	wildBase,
	variantBase
) {
	//////////////////////
	// Parameter Checks //
	//////////////////////

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

	//////////////////////
	//  Function Logic  //
	//////////////////////

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
		return {
			bestSnapbackBase: NUCLEOTIDE_COMPLEMENT[wildBase], // Snapback is matching wild
			bestDifference: wildSnapbackTmDiff,
		};
	} else {
		// Scenario B wins
		return {
			bestSnapbackBase: NUCLEOTIDE_COMPLEMENT[variantBase], // Snapback is matching variant
			bestDifference: variantSnapbackTmDiff,
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
	// 1) Validate Parameters
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
	baseUrl += `&decimalplaces=${SNAP_DECIMAL_PLACES}`;
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
 * Constructs or refines the snapback stem on the chosen strand so that
 * the melting temperature of the wild type allele gets as close to the passed
 * in desired snapback melting temperature for the wilt type
 *
 * @typedef {Object} StemMeltingTemp
 * @property {number} wildTm    - The snapbacks's melting temperature when extended with the wild allele
 * @property {number} variantTm - The snapbacks's melting temperature when matching the variant allele
 *
 * @typedef {Object} StemLoc
 * @property {number} start - The starting index (0-based) of the stem region in the sequence.
 * @property {number} end   - The ending index (0-based, inclusive) of the stem region in the sequence.
 *
 * @typedef {Object} CreateStemReturn
 * @property {StemLoc} stemLoc          - The finalized stem location in this strand context.
 * @property {StemMeltingTemp} meltingTemps - Melting temperatures for the wild and variant snapbacks
 * @property {number} snapbackBaseAtSNV - the base to use in the snapback (either complementary to the wild or variant base)
 *
 *
 * @param {string} targetStrandSeqSnapPrimerRefPoint - The DNA sequence (5'→3') for the chosen strand.
 * @param {Object} snvSiteSnapPrimerRefPoint         - SNV site object for this strand’s coordinates:
 *   - index: (number) The SNV's position in this strand.
 *   - variantBase: (string) The variant base (A/T/C/G).
 * @param {string} snapbackBaseAtSNV                 - The base (wild or variant) that the snapback is matching at the SNV.
 * @param {number} desiredSnapbackMeltTempWildType   - The target melting temperature (°C) for the wild-type snapback.
 *
 * @returns {CreateStemReturn} An object containing the finalized stem location and the calculated melting temperatures.
 */
async function createStem(
	targetStrandSeqSnapPrimerRefPoint,
	snvSiteSnapPrimerRefPoint,
	snapbackBaseAtSNV,
	desiredSnapbackMeltTempWildType
) {
	// Parameter Checking
	// Function logic
	// create initial stem again
}

/****************************************************************/
/*********************** Helper Functions ***********************/
/****************************************************************/

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
	////////////* Parameter checking *////////////////
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

	////////////* Logic of function */////////////////
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
	///////////// Parameter Checking ///////////////
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

	///////////// Core Logic ///////////////

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
 * Calculates the difference in melting temperature (Tm) for a given stem sequence
 * between a fully matched duplex and a mismatch duplex at a specified position.
 *
 * The mismatch base provided is assumed to be the base in the snapback that does
 * not match the snapback tail at the mismatch position — it could be either
 * the wild type or variant base, depending on which was not used in the tail.
 *
 * Assumptions:
 * - The snapback tail sequence is oriented 5' to 3'
 * - Mismatch is at a specific position and affects only one base pairing
 *
 * @param {string} currentStemSeqOnTail - The DNA stem sequence on the snapback tail (5'→3').
 * @param {number} mismatchPos - The position (0-based) in the sequence where the mismatch occurs.
 * @param {string} mismatchBase - The base in the snapback (A/T/C/G) that causes the mismatch.
 *
 * @returns {Promise<number>} - The absolute difference in Tm (°C) between the matched and mismatched stem.
 *
 * @throws {Error} - If inputs are invalid or Tm calculation fails.
 */
async function calculateStemTmDiff(
	currentStemSeqOnTail,
	mismatchPos,
	mismatchBase
) {
	//////////// Parameter Checking ////////////

	// Check sequence validity
	if (!isValidDNASequence(currentStemSeqOnTail)) {
		throw new Error(`Invalid DNA sequence: "${currentStemSeqOnTail}"`);
	}

	// Check mismatchPos
	if (
		typeof mismatchPos !== 'number' ||
		!Number.isInteger(mismatchPos) ||
		mismatchPos < 0 ||
		mismatchPos >= currentStemSeqOnTail.length
	) {
		throw new Error(
			`mismatchPos must be an integer in range [0, ${
				currentStemSeqOnTail.length - 1
			}]. Received: ${mismatchPos}`
		);
	}

	// Check mismatchBase
	if (
		typeof mismatchBase !== 'string' ||
		mismatchBase.length !== 1 ||
		!VALID_BASES.has(mismatchBase)
	) {
		throw new Error(
			`mismatchBase must be a single character from A, T, C, or G. Received: "${mismatchBase}"`
		);
	}

	//////////// Core Logic ////////////

	// 1. Calculate matched Tm
	const matchedTm = await getStemTm(currentStemSeqOnTail);

	// 2. Construct mismatch object
	const mismatch = {
		position: mismatchPos,
		type: mismatchBase,
	};

	// 3. Calculate mismatched Tm
	const mismatchedTm = await getStemTm(currentStemSeqOnTail, mismatch);

	// 4. Return absolute Tm difference
	return Math.abs(matchedTm - mismatchedTm);
}

/****************************************************************/
/******************** DNA Utility Functions *********************/
/****************************************************************/

/**
 * Checks if a given DNA sequence is valid. i.e. it is a string that only contains
 * the characters A, T, C, or G (case-insensitive).
 *
 * @param {string} seqStrand - The DNA sequence to validate.
 * @returns {boolean} - True passed in string is a valid DNA sequence
 */
function isValidDNASequence(seqStrand) {
	// Must be a string
	if (typeof seqStrand !== 'string') return false;

	// Must only contain uppercase A, T, C, or G
	// Its normal to use const here as its block scoped and base doesn't change in the block
	// for any iteration
	for (const base of seqStrand) {
		if (!VALID_BASES.has(base)) {
			return false;
		}
	}
	return true;
}

/**
 * Returns the DNA complement of a given nucleotide sequence.
 * If the sequence passed in is not a valid sequence it throws an error.
 *
 * @param {string} seqStrand - A string representing the DNA sequence (e.g., "ATCG").
 * @returns {string} - The complementary DNA sequence (e.g., "TAGC").
 */
function complementSequence(seqStrand) {
	if (!isValidDNASequence(seqStrand)) {
		throw new Error(`Invalid DNA sequence: ${seqStrand}`);
	}

	return seqStrand
		.toUpperCase()
		.split('') // Splits into an array of single character strings
		.map((base) => NUCLEOTIDE_COMPLEMENT[base])
		.join(''); // rejoins array into a string of the complement bases
}

/**
 * Returns the reverse complement of a DNA sequence.
 * Throws an error if the sequence is invalid.
 *
 * @param {string} seqStrand - A DNA sequence (e.g., "ATCG").
 * @returns {string} - The reverse complement (e.g., "CGAT").
 */
function reverseComplement(seqStrand) {
	const complementStrand = complementSequence(seqStrand);
	return complementStrand.split('').reverse().join('');
}

/**
 * Returns the reverse complement for a mismatch site.
 * That is, it returns the complement bases and corrects for the sequence's
 * orientation, assuming that the new sequence starts with the 5' end.
 *
 * Assumptions:
 * - Assumes the variant bases are all valid (elements of {"C", "G", "A", "T"}).
 *
 * Typedefs:
 * @typedef {Object} SNVSite
 * @property {number} index - The index in the sequence where the variant occurs.
 * @property {string} variantBase - A one character string representing the variant base.
 *
 *
 * @param {SNVSite} snvSite - An object representing the single nucleotide variant site.
 * @param {number} seqLen - The length of the target sequence.
 *
 * @returns {SNVSite} - An object representing the single nucleotide variant site for the reverse complement sequence.
 */

function revCompSNV(snvSite, seqLen) {
	////////* Checking arguements */////////
	// Ensure snvSite is an object with exactly two keys: 'index' and 'variantBase'
	if (
		typeof snvSite !== 'object' ||
		snvSite == null ||
		Array.isArray(snvSite) || // an array is also an object, but snvSite can't be an array
		Object.keys(snvSite).length !== 2 ||
		!('index' in snvSite) ||
		!('variantBase' in snvSite)
	) {
		throw new Error(
			'Invalid SNV object: must contain only { index: number, variantBase: string }'
		);
	}
	// Validate index: must be an integer in [0, seqLen - 1]
	if (
		typeof snvSite.index !== 'number' ||
		!Number.isInteger(snvSite.index) ||
		snvSite.index < 0 ||
		snvSite.index >= seqLen
	) {
		throw new Error(
			`Invalid SNV index: ${
				snvSite.index
			} (must be an integer from 0 to ${seqLen - 1}(seqLen -1))`
		);
	}
	// Validate variantBase: must be a valid base (A, T, C, G), case-sensitive
	if (
		typeof snvSite.variantBase !== 'string' ||
		snvSite.variantBase.length !== 1 ||
		!VALID_BASES.has(snvSite.variantBase)
	) {
		throw new Error(
			`Invalid variant base: "${snvSite.variantBase}" (must be one of 'A', 'T', 'C', or 'G')`
		);
	}

	// Since revComplement sequence starts with 5' end and both are indexed starting at 0
	const revCompIndex = seqLen - snvSite.index - 1;

	const revCompVariantBase = NUCLEOTIDE_COMPLEMENT[snvSite.variantBase];

	return {
		index: revCompIndex,
		variantBase: revCompVariantBase,
	};
}

/**
 * Reverses a DNA sequence.
 *
 * This function returns the input sequence in reverse order.
 *
 * Assumptions:
 * - The sequence consists only of valid DNA bases and is a string.
 *
 * @param {string} seqStrand - A string representing the DNA sequence (e.g., "ATCG").
 * @returns {string} - The reversed sequence (e.g., "GCTA").
 *
 * @throws {Error} - If the input is not a string or contains invalid characters.
 */
function reverseSequence(seqStrand) {
	// Validate input
	if (!isValidDNASequence(seqStrand)) {
		throw new Error(
			`Invalid DNA sequence: "${seqStrand}". Must only contain A, T, C, or G and be a string`
		);
	}

	// Reverse the sequence
	return seqStrand.split('').reverse().join('');
}

/****************************************************************/
/*********************** Export Function ************************/
/****************************************************************/

export {
	// Primary function
	createSnapback,

	// Secondary functions
	useTargetStrandsPrimerForComplement,
	evaluateSnapbackMatchingOptions,
	getStemTm,
	createStem,

	// Helper/logic functions
	snvTooCloseToPrimer,
	buildMismatchSequenceForAPI,
	parseTmFromResponse,
	calculateStemTmDiff,

	// DNA utility functions
	isValidDNASequence,
	complementSequence,
	reverseComplement,
	revCompSNV,
	reverseSequence,

	// Constants
	SNV_BASE_BUFFER,
	NUCLEOTIDE_COMPLEMENT,
};
