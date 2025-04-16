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
const SNV_BASE_BUFFER = 3; // The number of matched bases required on either end of a mismatched SNV
const INNER_LOOP_NUMBER_OF_STRONG_BASE_MISMATCHES_REQUIRED = 2;
const END_OF_STEM_NUMBER_OF_STRONG_BASE_MISMATCHES_REQUIRED = 2;
const MINIMUM_TARGET_SNAPBACK_MELTING_TEMP = 40;
const MAXIMUM_TARGET_SNAPBACK_MELTING_TEMP = 80;
const SNAP_DECIMAL_PLACES = 2;
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
	// Make sure the SNV is not too close to either end of the primer that a proper stem cannot be formed
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

	/******* Creating complementary strand's variables and initial, could-be snapback primer stems *********/
	/*** Getting complementary strand variable ***/
	// We must reverse the complement strand so that it starts with the 5' end
	const compTargetSeqStrand = reverseComplement(targetSeqStrand);
	const compSnvSite = revCompSNV(snvSite, targetSeqStrand.length);

	// Adds {SNV base buffer} matching neucleotides on each end of SNV so that mismatch is definitely included in stem
	const initStemLoc = {
		start: snvSite.index - SNV_BASE_BUFFER,
		end: snvSite.index + SNV_BASE_BUFFER,
	};
	// Doing the same for the complementary strand
	const compInitStemLoc = {
		start: compSnvSite.index - SNV_BASE_BUFFER,
		end: compSnvSite.index + SNV_BASE_BUFFER,
	};

	/* 
		Calculate the initial melting temperature differences for variants and choose the primer and base match (variant or wild)
		with the largest differnce as snapback 
	*/
	// Note: It is assumed that the melting temperature difference will not vary as the stem length is increased, in accordance
	//       with the nearest neighbor model of the Santa Lucia melting temperature equations
	const { useTargetStrand, snapbackBaseAtSNV } =
		await useTargetStrandsPrimerForComplement(
			targetSeqStrand,
			compTargetSeqStrand,
			initStemLoc,
			compInitStemLoc,
			snvSite,
			compSnvSite
		);

	// Assigning variables in terms of the primer strand to use as the snapback
	var targetStrandSeqSnapPrimerRefPoint;
	var snvSiteSnapPrimerRefPoint;
	var stemStartSnapPrimerRefPoint; // index for stem start in the primers reference point (tail is at target sequences 5' end with index 0)
	var stemEndSnapPrimerRefPoint;
	var allowedStemStart; //The stem can not go past this location or it will me complimentary to where the snapback primes to the target sequence
	var allowedStemEnd; // The stem can not go past this location or it will be complementary to the limiting primer
	var currLoopLen;

	if (useTargetStrand) {
		targetStrandSeqSnapPrimerRefPoint = targetSeqStrand;
		snvSiteSnapPrimerRefPoint = snvSite;
		stemStartSnapPrimerRefPoint = initStemLoc.start;
		stemEndSnapPrimerRefPoint = initStemLoc.end;
		allowedStemStart = primerLen;
		allowedStemEnd = targetSeqStrand.length - compPrimerLen - 1;
		// +2 accounts for the 2 base mismatch to prevent loop zipping
		currLoopLen = initStemLoc.start + 2;
	} else {
		targetStrandSeqSnapPrimerRefPoint = compTargetSeqStrand;
		snvSiteSnapPrimerRefPoint = compSnvSite;
		stemStartSnapPrimerRefPoint = compInitStemLoc.start;
		stemEndSnapPrimerRefPoint = compInitStemLoc.end;
		allowedStemStart = compPrimerLen;
		allowedStemEnd = compTargetSeqStrand.length - primerLen - 1;
		// +2 accounts for the 2 base mismatch to prevent loop zipping
		currLoopLen = compInitStemLoc.start + 2;
	}

	const objQuestionMark = createStem(
		targetStrandSeqSnapPrimerRefPoint,
		snvSiteSnapPrimerRefPoint,
		snapbackBaseAtSNV,
		stemStartSnapPrimerRefPoint,
		stemEndSnapPrimerRefPoint,
		allowedStemStart,
		allowedStemEnd,
		desiredSnapbackMeltTempWildType
	);
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
 * @param {string} compTargetSeqStrand - Full complementary strand (5'→3') [reverse complement of the target]
 * @param {Object} initStemLoc         - { start: number, end: number } for the target’s initial stem
 * @param {Object} compInitStemLoc     - { start: number, end: number } for the complementary’s initial stem
 * @param {Object} snvSite             - { index: number, variantBase: string }
 * @param {Object} compSnvSite         - { index: number, variantBase: string }
 *
 * @returns {Promise<{ useTargetStrand: boolean, snapbackBaseAtSNV: string }>}
 */
async function useTargetStrandsPrimerForComplement(
	targetSeqStrand,
	compTargetSeqStrand,
	initStemLoc,
	compInitStemLoc,
	snvSite,
	compSnvSite
) {
	// 1) Slice out the "init stem" region from each strand
	const targetInitStem = targetSeqStrand.slice(
		initStemLoc.start,
		initStemLoc.end
	);
	const compInitStem = compTargetSeqStrand.slice(
		compInitStemLoc.start,
		compInitStemLoc.end
	);

	// 2) Identify the wild-type base on each strand
	//    Because there's only one variant, it's the difference from the wild base.
	const targetWildBase = targetSeqStrand[snvSite.index];
	const compWildBase = compTargetSeqStrand[compSnvSite.index];

	// 3) SNV is at position 3 in these 7-base slices (SNV_BASE_BUFFER=3)
	const mismatchPos = 3;

	// 4) Evaluate Tm differences for the target strand
	const targetScenario = await evaluateSnapbackOptions(
		targetInitStem,
		mismatchPos,
		targetWildBase,
		snvSite.variantBase
	);

	// 5) Evaluate Tm differences for the complementary strand
	const compScenario = await evaluateSnapbackOptions(
		compInitStem,
		mismatchPos,
		compWildBase,
		compSnvSite.variantBase
	);

	// 6) Compare which scenario yields the bigger Tm difference
	if (targetScenario.bestDifference > compScenario.bestDifference) {
		return {
			useTargetStrand: true,
			snapbackBaseAtSNV: targetScenario.bestBase,
		};
	} else {
		return {
			useTargetStrand: false,
			snapbackBaseAtSNV: compScenario.bestBase,
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
 * @param {string} initStem     - The sub-sequence (7 bases) containing the SNV at position `mismatchPos`
 * @param {number} mismatchPos  - Usually 3 (SNV in center)
 * @param {string} wildBase     - The wild-type base at that SNV
 * @param {string} variantBase  - The single variant base
 *
 * @returns {Promise<{ bestBase: string, bestDifference: number }>}
 *   bestBase => the base (wild or variant) to use in the snapback
 *   bestDifference => the numeric Tm difference from the wild scenario
 */
async function evaluateSnapbackOptions(
	initStem,
	mismatchPos,
	wildBase,
	variantBase
) {
	// 1) Tm with wild base matched
	const wildMatchTm = await getStemTm(initStem);

	// 2) Tm with variant base matched

	// 2) Scenario A: Snapback matches the wild base => mismatch is the variant base
	// The base on the opposite strand for the mismatch will be the complement
	const mismatchObjA = {
		position: mismatchPos,
		type: NUCLEOTIDE_COMPLEMENT[variantBase],
	};
	const mismatchTmA = await getStemTm(initStem, mismatchObjA);
	const diffA = Math.abs(wildTm - mismatchTmA);

	// 3) Scenario B: Snapback matches the variant base => mismatch is the wild base
	const mismatchObjB = { position: mismatchPos, type: wildBase };
	const variantMatchTmB = await getStemTm(initStem, mismatchObjB);
	const diffB = Math.abs(wildTm - variantMatchTmB);

	// 4) Pick whichever scenario yields the larger difference
	if (diffA > diffB) {
		// Scenario A wins
		return {
			bestBase: wildBase, // Snapback is matching wild
			bestDifference: diffA,
		};
	} else {
		// Scenario B wins
		return {
			bestBase: variantBase, // Snapback is matching variant
			bestDifference: diffB,
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
	// 1) Validate input
	if (!seq || !isValidDNASequence(seq)) {
		throw new Error(`Invalid or empty sequence: "${seq}"`);
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
	let baseUrl = 'https://dna-utah.org/ths/cgi-bin/tmsnap.cgi';
	baseUrl += `?mg=${MG}`;
	baseUrl += `&mono=${MONO}`;
	baseUrl += `&seq=${seq}`;
	baseUrl += `&tparam=${T_PARAM}`;
	baseUrl += `&saltcalctype=${SALT_CALC_TYPE}`;
	baseUrl += `&otype=${O_TYPE}`;
	baseUrl += `&concentration=${CONC}`;
	baseUrl += `&limitingconc=${LIMITING_CONC}`;
	baseUrl += `&decimalplaces=${SNAP_DECIMAL_PLACES}`;
	// If a mismatch is passed in this will add the correct mismatch sequence
	if (mismatch) {
		baseUrl += `&mmseq=${mismatchSeq}`;
	}

	// 4) For local dev only. I think we just use the baseURL if this code is on server?
	const proxyUrl = `https://api.allorigins.win/get?url=${baseUrl}`;

	// 5) Fetch and parse response
	const response = await fetch(proxyUrl);
	if (!response.ok) {
		throw new Error(
			`Network error: ${response.status} - ${response.statusText}`
		);
	}

	const data = await response.json();
	if (!data || !data.contents) {
		throw new Error("Response missing 'contents'. Possibly a proxy error.");
	}

	const rawHtml = data.contents;

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
 * @typedef {Object} StemLoc
 * @property {number} start - The final start index of the stem on this strand.
 * @property {number} end   - The final end index of the stem on this strand.
 *
 * @typedef {Object} StemMeltingTemp
 * @property {number} wildTm    - The stem's melting temperature when extended with the wild allele
 * @property {number} variantTm - The stem's melting temperature when matching the variant allele
 *
 * @typedef {Object} CreateStemReturn
 * @property {StemLoc} stemLoc          - The finalized stem location in this strand context.
 * @property {StemMeltingTemp} meltingTemp - Melting temperatures for the wild and variant stems.
 *
 * Recursively constructs or refines the snapback stem on the chosen strand so that
 * its melting temperature meets the specified thresholds for both the wild and
 * variant bases at the SNV site.
 *
 * @param {string} targetStrandSeqSnapPrimerRefPoint - The DNA sequence (5'→3') for the chosen strand.
 * @param {Object} snvSiteSnapPrimerRefPoint         - SNV site object for this strand’s coordinates:
 *   - index: (number) The SNV's position in this strand.
 *   - variantBase: (string) The variant base (A/T/C/G).
 * @param {string} snapbackBaseAtSNV                 - The base (wild or variant) that the snapback is matching at the SNV.
 * @param {number} stemStartSnapPrimerRefPoint       - The initial start index of the stem in this strand.
 * @param {number} stemEndSnapPrimerRefPoint         - The initial end index of the stem in this strand.
 * @param {number} allowedStemStart                  - The furthest upstream index to which the stem can safely extend.
 * @param {number} allowedStemEnd                    - The furthest downstream index to which the stem can safely extend.
 * @param {number} desiredSnapbackMeltTempWildType   - The target melting temperature (°C) for the wild-type snapback.
 *
 * @returns {CreateStemReturn} An object containing the finalized stem location and the calculated melting temperatures.
 */
async function createStem(
	targetStrandSeqSnapPrimerRefPoint,
	snvSiteSnapPrimerRefPoint,
	snapbackBaseAtSNV,
	stemStartSnapPrimerRefPoint,
	stemEndSnapPrimerRefPoint,
	allowedStemStart,
	allowedStemEnd,
	desiredSnapbackMeltTempWildType
) {}

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
	// Must be at least {buffer} bases away from either primer
	const lowerBoundIndex = primerLen + buffer;
	const upperBoundIndex = seqLen - compPrimerLen - buffer - 1;

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
	///////////* Parameter Checking */////////////
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

	// Replace that position with the mismatch complement's in the original sequence, thisis what the API wants
	const arr = seq.split('');
	arr[mismatch.position] = complementBase;
	return arr.join('');
}

/**
 * Simple parser to extract the numeric Tm from the raw HTML:
 * e.g. <html><head></head><body><seq>...</seq><tm>47.27</tm><mmtm>37.54</mmtm></body></html>
 * ******I let chatGPT make this one. I'll go through it later and clean it up if required
 *
 * @param {string} rawHtml - The raw string returned by tmsnap.cgi
 * @param {boolean} [mismatch] - Optional specification denoting if we want to parse out the mismatched tm
 * @returns {number|null}  - The Tm if found, otherwise null
 */
function parseTmFromResponse(rawHtml, mismatch) {
	try {
		// Parse the text into a DOM
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

/****************************************************************/
/*********************** Export Function ************************/
/****************************************************************/

export {
	// Primary function
	createSnapback,

	// Secondary functions
	useTargetStrandsPrimerForComplement,
	evaluateSnapbackOptions,
	getStemTm,
	createStem,

	// Helper/logic functions
	snvTooCloseToPrimer,
	buildMismatchSequenceForAPI,
	parseTmFromResponse,

	// DNA utility functions
	isValidDNASequence,
	complementSequence,
	reverseComplement,
	revCompSNV,
};
