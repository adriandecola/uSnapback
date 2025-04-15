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
 *     - Passed in target sequence is valid and uppercase
 * - Assumes that the target sequence begins with the 5' end
 * - we want to minimize loop length for whatever reason??
 *
 *
 * Things to consider/add/make better:
 *
 *
 * @param {string} targetSeqStrand - A strand of the full DNA sequence to design the snapback primer for.
 * @param {number} primerLen - The length of the primer on the strand denoted by the target sequence
 * @param {number} compPrimerLen - The length of the complementary primer on the opposite strand.
 * @param {Object} snvSite - An object representing the single neucleotide variant site with:
 *     @property {number} index - The index in the sequence where the variant occurs (0-based).
 *     @property {string} variantBase - A one character string representing the variant base.
 * @param {number} minSnapbackMeltTemp - The minimum viable snapback melting temperature for any match or mismatch snapback
 * @param {number} desiredSnapbackMeltTemp = The desired snapback melting temperature for the matched snapback (no mismatches)
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
function createSnapback(
	targetSeqStrand,
	primerLen,
	compPrimerLen,
	snvSite,
	desiredSnapbackMeltTemp
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
	compTargetSeqStrand = reverseComplement(targetSeqStrand);
	compSnvSite = revCompSNV(snvSite, targetSeqStrand.length);

	// Adds {SNV base buffer} matching neucleotides on each end of SNV so that mismatch is definitely included in stem
	initStemLoc = {
		start: snvSite.index - SNV_BASE_BUFFER,
		end: snvSite.index + SNV_BASE_BUFFER,
	};
	// Doing the same for the complementary strand
	compInitStemLoc = {
		start: compSnvSite.index - SNV_BASE_BUFFER,
		end: compSnvSite.index + SNV_BASE_BUFFER,
	};

	/* 
		Calculate the initial melting temperature differences for variants and choose the primer and base match (variant or wild)
		with the largest differnce as snapback 
	*/
	// Note: It is assumed that the melting temperature difference will not vary as the stem length is increased, in accordance
	//       with the nearest neighbor model of the Santa Lucia melting temperature equations

	//
	// Future stem cannot not overlap with anywhere a primer would go
	/*****************only going to use one of these aassign to new variables */
	const allowedStemSeq = targetSeqStrand.slice(
		primerLen,
		targetSeqStrand.length - compPrimerLen
	);
	// Makes sure stem does not match with anywhere a primer would go
	compAllowedStemSeq = targetSeqStrand.slice(
		compPrimerLen,
		compTargetSeqStrand.length - primerLen
	);

	// As we do not want the loop length to go larger than
	const loopLen = initStemLoc.start + 2;

	snapbackPrimerTargetStrand = calculateStem(
		allowedStemSeq,
		initStemLoc,
		desiredSnapbackMeltTemp
	);
}

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
	// 1) Tm with wild base matched (i.e., no mismatch in initStem)
	const wildTm = await getStemTm(initStem);

	// 2) Scenario A: Snapback matches the wild base => mismatch is the variant base
	const mismatchObjA = { position: mismatchPos, type: variantBase };
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
 * This function creates a snapback primer given an input sequence that includes both the 5' end and 3' end
 * primers and the location and type of the single base variant of interest. The 'scorpion tail' of the snapback
 * primer is added to the 5' end. It includes a loop and then a tail that contains a bases complimentary to the
 * input sequence, allowing it to fold on itself (once duplicated) creating a hairpin structure.
 * The hybridized part includes a single base variant of interest. Snapback primers
 * are helpful because the melting temperature of the hairpin structure can vary more greatly for a single base
 * mismatch than it can for regular PCR product.
 *
 * This recursive function finds a snapback primer such that difference in melting temperatures for the wild
 * and variant type single base variant is maximized. For cases where the variant can be more than one base,
 * the mimimum temperature difference between any trait is maximized.
 *
 * Assumptions:
 * - There is only one mismatch site
 * - The stem (either end) should not be complementary to primer end
 *      - This means that the
 * - The hairpin loop must be at least 6 bases long; however this is satisfied as the primer should
 *      -> If the loop of the hairpin is too short, the model for melting temperatures is not accurate and it is hard
 *         for the hairpin structure to form
 *      -> however this assumption is satifies as the primer is always 15 base pairs long and the
 * - The stem's melting temperatuere is calcualted using the Santa Lucia method.
 * - Extension can occur on the hairpin's complement, on its 5' end. A 2 base pair mismatch after the stem can
 *   help avoid this
 * - The acceptable snapback melting temperature range is 50-75 degrees Celcius
 *
 * Process:
 * - The function starts by checking the termination conditions. It checks if the built up sequence
 *   we are testing for goes past the viable stem area  what if it goes over limiting primer end?????
 *   or if the stem is longer than **** base pairs. If this is the case the function returns the
 *   sequence with the largest temperature difference after adding a 2 base pair mismatch after the stem on the 5'
 *   end to avoid the snapbacks compliment from extending.
 * - If the function has no built up sequence yet, the function initiallizes with adding 3 complimentary
 *   neucleotices on each end of the varying neucleotide. This is only done if it has not been done before.
 *      ***is it done for both sides?
 * - The function then calculates the melting temperatures of the matched and mismatched sequence(s).
 * - If the melting temperatures of the matched and mismatched sequence(s) are inside greater than the mimimum
 *   snapback temperature of 50 degrees Celsius, then the function calculates the temperature difference
 *   (for snapbacks that are only testing for one possible base change) or the minimmum temperature
 *   difference between all different bases (for snapbacks that are testing for multiple possible
 *   base changes)
 * - If not then the snapback calls itself on the build-up sequence with added neucleotides on either end (as long as
 *   they do not intrude on either primer ends (don't wasn self annealing at those locations) to make the stem longer,
 *   so that the melting temperatures are in an acceptable range
 * - If this melting temperture difference is greater than any previously saved melting temperature
 *   difference, or if one has not been saved yet, the melting temperature difference is saved.
 * - The function calls its self of the builtup sequence with all the possible addition types of
 *
 *
 *
 *
 * Notes to make this equation more effecient:
 * - I could represent DNA sequences as an array of 2 bit encoded neucleotides. This could save some memory
 * - I could hash and cache melting temperatures of stems? ACTUALLY they shouldnt repeat much unless theres convenient patterns in the sequence so nevermind? or try it and see how many times a cached item is used
 * - I could save the melting temperatur so far and implement santa lucia method here and simply add entropys
 *   and enthalpies of added neucleotides.
 * - have a clean function to call a more simple recursive function to take out some cleaning each call
 * - make my own data type to represent dna sequences with their own methods/properties like .complement or .complement()
 *
 *
 * ** it doesnt check if it overlaps with primers or anything, but seems very rare, can force it not to if loop should be extended (add 2 pair mismatch?)
 *
 *
 *
 * @param {string} targetSeq - The DNA sequence to design the primer (snapback) for. Assume limiting
 *                             primer is on complementary sequence.
 *
 * @param {number} primerLen - The length of the primer (excluding snapback). Assume limiting
 *                             primer is on complementary sequence.
 * @param {number} compPrimerLen - The length of the complementary primer.
 * @param {number} minLoopLen - The minimum loop length so far, assuming that the entire tail is part
 *                              of the stem. Includes a 2 base mismatch at start of loop towars 5' end
 *                              to make sure loop is not self complimentary and closes itself
 * @param {string} builtUpSeq - The DNA sequence that has been built up for testing.
 * @param {string} allowedStemSeq - The section of the DNA sequence that can be used as part of the stem
 * @param {number} minimumLoopLen -> this is the primer length on that end and some extra
 * @returns {string} - The full snapback primer sequence (5'-tail+primer-3').
 * @property {}
 * shoult return the loop size, and entire sequence, and stem size?
 * ** should also return the  matched temperature?? and the mismatched temperature
 ***dont need 2 base pair mismatch every time
 */
function createStem() {}

/****************************************************************/
/*********************** Helper Functions ***********************/
/****************************************************************/

/**
 * Checks if a given DNA sequence is valid. i.e. it is a string that only contains
 * the characters A, T, C, or G (case-insensitive).
 *
 * @param {string} seqStrand - The DNA sequence to validate.
 * @returns {boolean} - True passed in string is a valid DNA sequence
 */
function isValidDNASequence(seqStrand) {
	/* 
    'const' is used as 'base' should not change within an iteration. It still changes
    every successive iterations though. 
    */
	for (const base of seqStrand.toUpperCase()) {
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
	// Since revComplement sequence starts with 5' end and both are indexed starting at 0
	revCompIndex = seqLen - snvSite.index - 1;

	revCompVariantBase = NUCLEOTIDE_COMPLEMENT[snvSite.variantBase];

	return {
		index: revCompIndex,
		variantBase: revCompVariantBase,
	};
}

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
 */
function snvTooCloseToPrimer(snvIndex, primerLen, compPrimerLen, seqLen) {
	// Must be at least {buffer} bases away from either primer
	const lowerBoundIndex = primerLen + buffer;
	const upperBoundIndex = seqLen - compPrimerLen - buffer - 1;

	return snvIndex < lowerBoundIndex || snvIndex > upperBoundIndex;
}

/**
 * Calculates the melting temperature of a snapback stem via server API (dna-utah.org).
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
	if (mismatchSeq) {
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
	const tmValue = parseTmFromResponse(rawHtml);

	if (tmValue === null) {
		throw new Error(
			'No <tm> element found or invalid numeric value in server response.'
		);
	}

	// 6) Return the numeric Tm
	return tmValue;
}

/**
 * Simple parser to extract the numeric Tm from the raw HTML:
 * e.g. <html><head></head><body><seq>...</seq><tm>47.27</tm><mmtm>37.54</mmtm></body></html>
 * ******I let chatGPT make this one. I'll go through it later and clean it up if required
 *
 * @param {string} rawHtml - The raw string returned by tmsnap.cgi
 * @returns {number|null}  - The Tm if found, otherwise null
 */
function parseTmFromResponse(rawHtml) {
	try {
		// Parse the text into a DOM
		const parser = new DOMParser();
		const doc = parser.parseFromString(rawHtml, 'text/html');

		// Look for a <tm> element
		const tmElement = doc.querySelector('tm');
		if (!tmElement) {
			return null;
		}

		// Convert the text inside <tm> to a float
		const tmValue = parseFloat(tmElement.textContent.trim());
		return isNaN(tmValue) ? null : tmValue;
	} catch (err) {
		console.error('parseTmFromResponse error:', err);
		return null;
	}
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
 * @property {string} type     - The base you want on the opposite strand (e.g. "G").
 *
 *
 * @param {string} seq - Original matched sequence (5'→3').
 * @param {Mismatch} mismatch - The mismatch specification.
 *
 * @throws {Error} If the sequence or mismatch is invalid.
 * @returns {string} The mmseq string for the Tm service.
 */
function buildMismatchSequenceForAPI(seq, mismatch) {
	// Ensure the main sequence is valid. Throws an error if invalid or empty.
	isValidDNASequence(seq);

	// Validate mismatch object
	if (
		!mismatch ||
		mismatch.position == null ||
		typeof mismatch.type !== 'string'
	) {
		throw new Error(`Invalid mismatch object: ${JSON.stringify(mismatch)}`);
	}

	// Make sure the mismatch.type is a single valid base (A/T/C/G).
	// We'll treat mismatch.type as a minimal "sequence" to check:
	isValidDNASequence(mismatch.type);

	// Ensure position is within seq bounds
	if (mismatch.position < 0 || mismatch.position >= seq.length) {
		throw new Error(
			`Mismatch position ${mismatch.position} is out of range for sequence length ${seq.length}.`
		);
	}

	// The Tm service wants the difference in `mmseq` to be the complement
	// of the mismatch base. (e.g., mismatch.type='G' -> insert 'C').
	const complementBase = NUCLEOTIDE_COMPLEMENT[mismatch.type.toUpperCase()];
	// If mismatch.type wasn't in {A,T,C,G}, isValidDNASequence would have already thrown,
	// but we check again out of caution.
	if (!complementBase) {
		throw new Error(
			`Mismatch type "${mismatch.type}" has no complement. Should never happen.`
		);
	}

	// Replace the character at 'position' with that complement
	const arr = seq.split('');
	arr[mismatch.position] = complementBase;

	return arr.join('');
}
