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
const VALID_BASES = new Set(['A', 'T', 'C', 'G']);
//The number of matched bases required on either end of a mismatched SNV for snapback primers
const SNV_BASE_BUFFER = 3; 
const MINIMUM_SNAPBACK_LOOP_LENGTH = 6;
const MINIMUM_TARGET_SNAPBACK_MELTING_TEMP = 40;
const MAXIMUM_TARGET_SNAPBACK_MELTING_TEMP = 80;
// Min and max snapback melting temperature sthe users should be able to input shoudl be like 40-90?

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

/**
 * Creates a snapback primer sequence by identifying suitable stem regions
 * in the target DNA sequence, accounting for primer lengths and a single
 * mismatch site. This function initializes the process by preparing the
 * allowed stem sequence and then calls a recursive helper function to find
 * the best snapback configuration.
 *
 * Assumptions:
 * - Assumes passed in values are correct
 *     - Passed in target sequence is valid and uppercase
 * - Assumes that the target sequence begins with the 5' end 
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
 *     @property {Array<string>} variantBases - An array of one or more of the possible variant bases
 * @param {number} minSnapbackMeltTemp - The minimum viable snapback melting temperature for any match or mismatch snapback
 * @param {number} desiredSnapbackMeltTemp = The desired snapback melting temperature for the matched snapback (no mismatches)
 *
 * @returns {Object} - An object representing the formed snapback primer
 *     @property {string} snapbackSeq - The final snapback primer sequence (5'-tail + primer-3').
 *     @property {boolean} isOnTargetPrimer - A boolean representing if the primer to add a tail on is the primer denoted
 *                                            by the target sequence (True) or by its complement (False)
 *     @property {Array<Object>} snapbackMeltingTemps - Melting temperature info for match and mismatch variants.
 *         @property {string} snapbackMeltingTemps[].base - The base used at the mismatch site (A, T, C, or G).
 *         @property {boolean} snapbackMeltingTemps[].isMatch - True if this base is the correct one (no mismatch in snapback).
 *         @property {number} snapbackMeltingTemps[].meltingTemp - The melting temperature for this base configuration.
 *
 */
function createSnapback(
	targetSeqStrand,
	primerLen,
	compPrimerLen,
	snvSite,
	minSnapbackMeltTemp,
	desiredSnapbackMeltTemp
) {
    // Makes sure stem does not match with anywhere a primer would go
	const allowedStemSeq = targetSeqStrand.slice(
		primerLen,
		targetSeqStrand.length - compPrimerLen
	);

    // Adds 3 matching neucleotides on each end of SNV so that mismatch is definitely included in stem
	currStemLoc = {
		start: snvSite.index - 3,
		end: snvSite.index + 3,
	};

	// +2 for the two base mismatches in the loop to avoid loop closing
	// (assumes the tail with the 2 base pairs added is part of the loop)
	// more loop can be added after the fact
	const minLoopLen = currStemLoc.start + 2;

	snapbackPrimerTargetStrand = calculateStem(
		allowedStemSeq,
		currStemLoc,
		minLoopLen,
		minSnapbackMeltTemp,
		desiredSnapbackMeltTemp
	);

    /******** Doing the same for the complement strand ********/

    // We must reverse the complement strand so that it starts with the 5' end 
    compTargetSeqStrand = reverseComplement(targetSeqStrand);
    
    // Mismatch Site 

    // Makes sure stem does not match with anywhere a primer would go
    compAllowedStemSeq = targetSeqStrand.slice(
        compPrimerLen,
        compTargetSeqStrand.length - primerLen
    )

    compCurrStemLoc = {
        start: 
    }







	snapbackPrimerCompStrand = calculateStem(
		,

	);

	// First build up the 3 padding on each end for the stem

	//
}





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
 * Returns the reverse compelment for a mismatch site. 
 * That is it returns the complement bases and correct for the sequence's
 * complement assuming that the new sequence starts with the 5' end
 * 
 * Assumptions:
 * - Assumes the variant bases are all valid (element of {"C", "G", "A", "T"})
 * 
 * 
 * @param {Object} snvSite - An object representing the single neucleotide variant site with:
 *     @property {number} index - The index in the sequence where the variant occurs.
 *     @property {Array<string>} variantBases - An array of one or more of the possible variant bases
 * @param {number} seqLen - The length of the target sequence
 * 
 * @returns {Object} - An object representing the single neucleotide variant site for the reverse complement sequence with:
 *     @property {number} index - The index in the sequence where the variant occurs.
 *     @property {Array<string>} variantBases - An array of one or more of the possible variant bases
 */
function revCompSNV(snvSite, seqLen) {
    // Since revComplement sequence starts with 5' end and both are indexed starting at 0
    revCompIndex = seqLen - snvSite.index - 1;


    revCompVariantBases = snvSite.variantBases.map(base => NUCLEOTIDE_COMPLEMENT[base]);

    return {
        index:revCompIndex,
        variantBases: revCompVariantBases
    }
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