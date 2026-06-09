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
const SNV_BASE_BUFFER = 3; // The number of matched bases required on either end of a mismatched SNV
const INNER_LOOP_NUMBER_OF_STRONG_BASE_MISMATCHES_REQUIRED = 1;
const END_OF_STEM_NUMBER_OF_STRONG_BASE_MISMATCHES_REQUIRED = 1;
const MINIMUM_TARGET_SNAPBACK_MELTING_TEMP = 40;
const MAXIMUM_TARGET_SNAPBACK_MELTING_TEMP = 80;
const MAX_AMPLICON_LEN = 1000;
const MIN_LOOP_LEN = 6;
const MIN_PRIMER_LEN = 12;
const TM_DECIMAL_PLACES = 2;
// Chemisty parameters
const MG = 2.2;
const MONO = 20.0;
const T_PARAM = 'SantaLuciaHicks';
const SALT_CALC_TYPE = 'bpdenominator';
const O_TYPE = 'oligo';
const PRIMER_O_TYPE = 'primer';
const CONC = 0.5;
const LIMITING_CONC = 0.5;
// This gets replaced by build.js
const API_URL = __API_URL__;
const PROXY_URL = __PROXY_URL__;
const USE_PROXY = __USE_PROXY__; // true  |  false  (a real Boolean)
const USE_TOKEN = __USE_TOKEN__; // true | false (real Boolean)
const API_TOKEN = __API_TOKEN__; // string or "" (build-time)

// Optional Santa Lucia/Rochester Tm data and calculators live in optionalTmMethods.js.
// Set this to false to run only the primary Wittwer method in the app.
const ENABLE_OPTIONAL_TM_METHODS = true;

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
 *        snapback-tail • optional inner-loop mismatch • stem (with chosen SNV base) • primer.
 *  4. Constructing descriptive objects for the unextended and extended products:
 *        - Unextended: a 5'→3' breakdown of the snapback primer itself.
 *        - Extended: segment strings taken directly from seq using the stem position as the guide:
 *            stuffBetween (includes the forward primer and any additional bases up to but not including the optional inner-loop block) •
 *            optional 5' inner-loop mismatch (left of the stem) •
 *            threePrimeStem (the stem interval) •
 *            threePrimerLimSnapExtMismatches (right of the stem) •
 *            threePrimerRestOfAmplicon (everything after).
 *  5. Annotating the SNV location in both extended descriptors:
 *        - snvOnThreePrimeStem.indexInThreePrimeStem is 0-based and relative to descriptiveExtendedSnapback.threePrimeStem.
 *        - snvOnFivePrimeStem.indexInFivePrimeStem is 0-based and relative to the 5' stem on that object
 *          (reverse-ordered relative to seq for the snapback side; derived by reversal on the limiting side).
 *        - For the snapback object, tailBaseAtSNV and matchesWild/matchesVariant indicate which allele the tail complements.
 *  6. Returning the snapback sequence, the limiting primer, Tm values and ΔTms, and the descriptive objects.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Assumptions
 * ──────────────────────────────────────────────────────────────────────────
 * - targetSeqStrand is a valid uppercase DNA string given 5'→3'.
 * - primerLen and compPrimerLen ≥ {MIN_PRIMER_LEN}.
 * - The SNV is ≥ {SNV_BASE_BUFFER} bases away from both primers.
 * - Stem Tm values are computed with the Santa Lucia nearest-neighbour model
 *   via dna-utah.org.
 * - We want to miminize the loop length to the primer on which the snapback tail is on.
 *   If the 5′ base of that primer complements the base immediately left of the stem,
 *   one strong mismatch is inserted at the 3′ end of the snapback tail so the loop
 *   does not zip.
 * - We want to keep the SNV in the middle of the stem or as close to centered as we
 *   can, as we build the snapback stem
 * - Extension can occur on the snapback's complement, on the side of the stem does not contain
 *   the loop, if the polymerase can easily attach to that location. A 2 base pair, strong mismatch,
 *   is added after the 5' end of the stem (on the snapback primer) help avoid this extension
 *   on its complement snapback
 * - All sequences are expressed 5'→3' in the frame of the primer that receives the tail.
 * - The optional inner-loop strong mismatch is immediately 5' of the stem (left of the stem).
 * - The strong end-of-stem mismatch used to block extension on the complementary snapback
 *   occur immediately 3' of the stem (right of the stem on this strand).
 * - For the extended descriptive object, stuffBetween includes the forward primer and all
 *   subsequent bases up to (but not including) the first inner-loop mismatch base.
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
 * @typedef {Object} SnapbackMeltingTempDiffs
 * @property {{matchWild:number, matchVariant:number}} onForwardPrimer  ΔTms (°C) if tail is on the forward primer
 * @property {{matchWild:number, matchVariant:number}} onReversePrimer  ΔTms (°C) if tail is on the reverse primer
 *
 * @typedef {Object} SnapbackMeltingTm
 * @property {number} 				wildTm     			Calculated Tm (°C) of the snapback on the wilt type allele
 * @property {number} 				variantTm  			Calculated Tm (°C) of the snapback on the variant type allele
 *
 * @typedef {Object} SnapbackPrimerResult
 * @property {string}						snapbackSeq			Entire snapback primer written 5' → 3'
 * @property {DescriptiveUnExtendedSnapbackPrimer} descriptiveUnExtendedSnapbackPrimer  Segment breakdown of the unextended snapback primer.
 * @property {DescriptiveExtendedSnapback}         descriptiveExtendedSnapback          Extended product segments and SNV indices (canonical threePrimeStem index).
 * @property {DescriptiveExtendedLimitingSnapback} descriptiveExendedLimSnapback        Limiting extended product segments and SNV indices.
 *                                                   			(tail → primer).
 * @property {string}						limitingPrimerSeq	The limiting primmer written 5' → 3'
 * @property {boolean}						tailOnForwardPrimer	true if tail is appended to the forward primer, i.e. the
 * 																primer represented by `targetSeqStrand`; false if it is
 * 																appended to the reverse primer.
 * @property {boolean}						matchesWild			true if the snapback base at the SNV matches on is tail
 *                                			       				the wild-type allele
 * @property {SnapbackMeltingTm}			snapbackMeltingTms	Object holding wild/variant snapback Tm values.
 * @property {SnapbackMeltingTempDiffs}		meltingTempDiffs 	Wild/variant ΔTm values for both primer orientations.
 *
 *
 * @typedef {Object} DescriptiveUnExtendedSnapbackPrimer
 * @property {string} fivePrimerLimSnapExtMismatches  Strong mismatch placed 3' of the stem on the snapback’s complement side (5' segment in snapback).
 * @property {string} fivePrimeStem                    Reverse-complement of seq[stem.start..stem.end], with the chosen tail base at the SNV.
 * @property {string} fivePrimeInnerLoopMismatches     Optional strong mismatch immediately 5' of the stem.
 * @property {string} forwardPrimer                    The forward primer sequence (seq.slice(0, primerLen)).
 *
 * @typedef {Object} SNVOnThreePrimeStem
 * @property {number} indexInThreePrimeStem  0-based index of the SNV relative to threePrimeStem.
 * @property {string} wildBase               Base in seq at the SNV index (wild-type).
 * @property {string} variantBase            Variant base at the SNV site.
 *
 * @typedef {Object} SNVOnFivePrimeStem
 * @property {number} indexInFivePrimeStem  0-based index of the SNV relative to fivePrimeStem for that object.
 * @property {string} tailBaseAtSNV         Base used on the snapback tail at the SNV (empty for limiting object).
 * @property {boolean} matchesWild          True if tailBaseAtSNV complements wildBase.
 * @property {boolean} matchesVariant       True if tailBaseAtSNV complements variantBase.
 * @property {string} compWildBase          Complement of wildBase.
 * @property {string} compVariantBase       Complement of variantBase.
 *
 * @typedef {Object} DescriptiveExtendedSnapback
 * @property {string} fivePrimerLimSnapExtMismatches   Strong mismatch placed to prevent extension on the complementary snapback.
 * @property {string} fivePrimeStem                    Snapback 5' stem segment (reverse-complement orientation).
 * @property {string} fivePrimeInnerLoopMismatches     Optional inner-loop strong mismatch immediately 5' of the stem.
 * @property {string} stuffBetween                     Sequence before the optional inner-loop block.
 * @property {string} threePrimeInnerLoopMismatches    Optional seq slice of the inner-loop block directly left of the stem.
 * @property {string} threePrimeStem                   seq.slice(stem.start, stem.end+1).
 * @property {string} threePrimerLimSnapExtMismatches  seq slice immediately right of the stem containing the extension-block mismatch.
 * @property {string} threePrimerRestOfAmplicon        Remainder of seq to the 3' end after the right-side mismatch block.
 * @property {SNVOnThreePrimeStem} snvOnThreePrimeStem SNV info indexed to threePrimeStem.
 * @property {SNVOnFivePrimeStem}  snvOnFivePrimeStem  SNV info indexed to fivePrimeStem (snapback side).
 *
 * @typedef {Object} DescriptiveExtendedLimitingSnapback
 * @property {string} threePrimerLimSnapExtMismatches
 * @property {string} threePrimeStem
 * @property {string} threePrimeInnerLoopMismatches
 * @property {string} stuffBetween
 * @property {string} fivePrimeInnerLoopMismatches
 * @property {string} fivePrimeStem
 * @property {string} fivePrimerLimSnapExtMismatches
 * @property {string} fivePrimerRestOfAmplicon
 * @property {SNVOnThreePrimeStem} snvOnThreePrimeStem SNV info using the same canonical index as the main object.
 * @property {SNVOnFivePrimeStem}  snvOnFivePrimeStem  SNV info indexed to fivePrimeStem on the limiting object
 *                                                     (tailBaseAtSNV is empty; indices derived by reversal).
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Parameters, Returns, and Errors
 * ──────────────────────────────────────────────────────────────────────────
 * @param {string}				targetSeqStrand					A strand of the full DNA sequence to design the snapback primer for.
 * 																(5'->3')
 * @param {number}				primerLen						The length of the forward primer (has the same bases as the beginning of
 * 																`targetSeqStrand`)
 * @param {number}				compPrimerLen					The length of the reverse primer
 * @param {SNVSite}				snvSite							An object representing the single nucleotide variant site
 * @param {number}				targetSnapMeltTemp				The desired snapback melting temperature for the wild type allele
 *
 * @returns {Promise<SnapbackPrimerResult>} 	 				Final snapback and limiting primers, snapback Tms and ΔTms,
 *                                           					and descriptive objects for unextended/extended products with SNV indices.
 *
 * @throws {Error} 												If any input is invalid, the SNV is too close to a primer,
 *             													or an acceptable stem cannot be constructed.
 */
async function createSnapback(
	targetSeqStrand,
	primerLen,
	compPrimerLen,
	snvSite,
	targetSnapMeltTemp,
) {
	//──────────────────────────────────────────────────────────────────────────//
	// Parameter Checking                                                      //
	//──────────────────────────────────────────────────────────────────────────//

	// 1. Validate the target sequence
	if (!isValidDNASequence(targetSeqStrand)) {
		throw new Error(
			`Invalid DNA sequence: "${targetSeqStrand}". Must be non-empty and contain only A, T, C, or G.`,
		);
	}

	// 2. Validate the desired snapback melting temperature
	if (
		typeof targetSnapMeltTemp !== 'number' ||
		!Number.isFinite(targetSnapMeltTemp) ||
		targetSnapMeltTemp <= 0
	) {
		throw new Error(
			'targetSnapMeltTemp must be a positive, finite number.',
		);
	}

	// 3. Validate primer lengths
	for (const [name, len] of [
		['primerLen', primerLen],
		['compPrimerLen', compPrimerLen],
	]) {
		if (typeof len !== 'number' || !Number.isInteger(len) || len < 0) {
			throw new Error(`${name} must be a non-negative integer.`);
		}
		if (len < MIN_PRIMER_LEN) {
			throw new Error(
				`${name} must be at least ${MIN_PRIMER_LEN} bases.`,
			);
		}
	}

	// 3a. Ensure primers can actually fit on the sequence
	if (primerLen + compPrimerLen >= targetSeqStrand.length) {
		throw new Error(
			`primerLen (${primerLen}) + compPrimerLen (${compPrimerLen}) ` +
				`cannot equal or exceed sequence length (${targetSeqStrand.length}).`,
		);
	}

	// 4. Validate the SNV object
	if (!isValidSNVObject(snvSite)) {
		throw new Error(
			`Invalid snvSite: ${JSON.stringify(
				snvSite,
			)}. Expected { index: number, variantBase: "A"|"T"|"C"|"G" }.`,
		);
	}
	if (snvSite.index >= targetSeqStrand.length) {
		throw new Error(
			`snvSite.index (${snvSite.index}) exceeds sequence length ${targetSeqStrand.length}.`,
		);
	}

	// 5. Ensure the SNV is sufficiently distant from both primers
	if (
		snvTooCloseToPrimer(
			snvSite.index,
			primerLen,
			compPrimerLen,
			targetSeqStrand.length,
		)
	) {
		throw new Error(
			`SNV at index ${snvSite.index} is too close to a primer; ` +
				`it must be at least ${SNV_BASE_BUFFER} bases away from both primer binding regions.`,
		);
	}

	// 6. Ensure the amplicon length does not exceed MAX_AMPLICON_LEN
	if (targetSeqStrand.length > MAX_AMPLICON_LEN) {
		throw new Error(
			`Amplicon length (${targetSeqStrand.length}) exceeds maximum allowed (${MAX_AMPLICON_LEN}).`,
		);
	}

	//──────────────────────────────────────────────────────────────────────────//
	//								Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1) Calculate the initial melting temperature differences for variants and choose the primer and base match (variant or wild),
	//	  with the largest difference, as the snapback primer.
	//
	// 	  Note: It is assumed that the melting temperature difference will scale(decrease) stem length is increased, but that the
	//	  snapback choice with the greatest temperature difference will remain the same bestSnapbackTailBaseAtSNV
	const {
		tailOnForwardPrimer,
		bestSnapbackTailBaseAtSNV: snapbackTailBaseAtSNV,
		snapbackTailMatchesWild: matchesWild,
	} = await useForwardPrimer(targetSeqStrand, snvSite);

	// 2) Assigning variables in terms of the primer strand to use as the snapback
	//	  It also creates variables for the information in the reverse complement frame of reference.
	//	  This is used later for calculating the melting temperature differences for
	var targetStrandSeqSnapPrimerRefPoint;
	var snvSiteSnapPrimerRefPoint;
	var primerLensSnapPrimerRefPoint;
	if (tailOnForwardPrimer) {
		targetStrandSeqSnapPrimerRefPoint = targetSeqStrand;
		snvSiteSnapPrimerRefPoint = snvSite;
		primerLensSnapPrimerRefPoint = {
			primerLen: primerLen,
			compPrimerLen: compPrimerLen,
		};
	} else {
		targetStrandSeqSnapPrimerRefPoint = reverseComplement(targetSeqStrand);
		snvSiteSnapPrimerRefPoint = revCompSNV(snvSite, targetSeqStrand.length);
		primerLensSnapPrimerRefPoint = {
			primerLen: compPrimerLen,
			compPrimerLen: primerLen,
		};
	}

	console.log('CREATING STEM_______________________');
	// 3) Calculating the stem in the snapback primers reference point (inculsive on both ends)
	const { bestStemLoc, meltingTemps } = await createStem(
		targetStrandSeqSnapPrimerRefPoint,
		snvSiteSnapPrimerRefPoint,
		primerLensSnapPrimerRefPoint,
		snapbackTailBaseAtSNV,
		matchesWild,
		targetSnapMeltTemp,
	);
	console.log('DONE CREATING STEM_______________________');

	// 4) Create a final snapback primer (in its reference point)
	const {
		snapback,
		descriptiveUnExtendedSnapbackPrimer,
		descriptiveExtendedSnapback,
		descriptiveExendedLimSnapback,
	} = buildSnapbackAndFinalProducts(
		targetStrandSeqSnapPrimerRefPoint,
		snvSiteSnapPrimerRefPoint,
		primerLensSnapPrimerRefPoint,
		bestStemLoc,
		snapbackTailBaseAtSNV,
	);

	const { snapbackTmRochester, snapbackTmSantaLucia } =
		await calculateOptionalSnapbackTms(descriptiveExtendedSnapback);

	// 5) Calculate melting temperature differences if we kept the same
	//	  stem location but changed the primer for which we attach the snapback
	//	  tail to or change which allele we match the nucleotide on the tail
	//	  to
	const meltingTempDiffs = await calculateMeltingTempDifferences(
		targetStrandSeqSnapPrimerRefPoint,
		snvSiteSnapPrimerRefPoint,
		bestStemLoc,
		tailOnForwardPrimer,
	);

	// 6) Return the results
	return {
		snapbackSeq: snapback,
		limitingPrimerSeq: reverseComplement(
			targetStrandSeqSnapPrimerRefPoint.slice(
				targetStrandSeqSnapPrimerRefPoint.length -
					primerLensSnapPrimerRefPoint.compPrimerLen,
			),
		),
		tailOnForwardPrimer: tailOnForwardPrimer,
		matchesWild: matchesWild,
		snapbackMeltingTms: meltingTemps,
		meltingTempDiffs: meltingTempDiffs,

		descriptiveUnExtendedSnapbackPrimer,
		descriptiveExtendedSnapback,
		descriptiveExendedLimSnapback,

		// Rochester hairpin model across loop + terminal mismatches + stem NN with/without SNV
		// This is the full object returned by calculateSnapbackTmRochester
		snapbackTmRochester,

		// SantaLucia hairpin model across loop + terminal mismatches + stem NN with/without SNV
		// This is the full object returned by calculateSnapbackTmSantaLucia
		snapbackTmSantaLucia,
	};
}

async function calculateOptionalSnapbackTms(descriptiveExtendedSnapback) {
	if (!ENABLE_OPTIONAL_TM_METHODS) {
		return {
			snapbackTmRochester: null,
			snapbackTmSantaLucia: null,
		};
	}

	const [snapbackTmRochester, snapbackTmSantaLucia] = await Promise.all([
		calculateSnapbackTmRochester(descriptiveExtendedSnapback),
		calculateSnapbackTmSantaLucia(descriptiveExtendedSnapback),
	]);

	return { snapbackTmRochester, snapbackTmSantaLucia };
}

async function calculateSnapbackTmRochester(descriptiveExtendedSnapback) {
	const optionalTmMethods = await import('./optionalTmMethods.js');
	return optionalTmMethods.calculateSnapbackTmRochester(
		descriptiveExtendedSnapback,
		getOptionalTmMethodOptions(),
	);
}

async function calculateSnapbackTmSantaLucia(descriptiveExtendedSnapback) {
	const optionalTmMethods = await import('./optionalTmMethods.js');
	return optionalTmMethods.calculateSnapbackTmSantaLucia(
		descriptiveExtendedSnapback,
		getOptionalTmMethodOptions(),
	);
}

function getOptionalTmMethodOptions() {
	return { getThermoParams, conc: CONC, limitingConc: LIMITING_CONC };
}

/*****************************************************************************************/
/********************************** Secondary Functions **********************************/
/*****************************************************************************************/

/**
 * Calculates melting temperature differences for alternative snapback primer
 * configurations, using the same stem location identified during snapback
 * construction. Specifically, it evaluates how the snapback Tm values differ
 * when:
 *   - The snapback tail is attached to the opposite primer strand
 *   - The snapback base at the SNV is matched to either the wild-type or variant allele
 *
 * This allows validation and comparison of snapback performance under different
 * configurations, ensuring robustness of the chosen snapback design.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Assumptions
 * ──────────────────────────────────────────────────────────────────────────
 * - Inputs have already been validated by `createSnapback`:
 *   - `targetStrandSeqSnapPrimerRefPoint` is a valid uppercase DNA sequence (A/T/C/G only).
 *   - `snvSiteSnapPrimerRefPoint` is a valid SNV object with `index` and `variantBase`.
 *     correspond the variant index and nucleotide.
 *   - `bestStemLoc` makes sense and does not go out of bounds of the sequence. It contains
 * 		only start and end keys that are valid integers
 *   - `tailOnForwardPrimer` is a boolean
 * - The mismatch objects used in calculating Tm are derived from valid
 *   wild/variant bases and aligned to the correct positions in their stem
 *   sequences.
 * - The Santa Lucia nearest-neighbour model (via `calculateSnapbackTmWittwer`) is
 *   used to compute Tm values.
 *
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Type Definitions
 * ──────────────────────────────────────────────────────────────────────────
 * @typedef {Object} MeltingTempDiffs
 * @property {Object} onForwardPrimer
 * @property {number|null} onForwardPrimer.matchWild    Absolute Tm difference (between wild and variant allele)
 * 														when tail is on forward primer, matching wild allele
 * @property {number|null} onForwardPrimer.matchVariant	Absolute Tm difference (between wild and variant allele)
 * 														when tail is on forward primer, matching variant allele
 * @property {Object} onReversePrimer
 * @property {number|null} onReversePrimer.matchWild    Absolute Tm difference (between wild and variant allele)
 * 														when tail is on reverse primer, matching wild allele
 * @property {number|null} onReversePrimer.matchVariant Absolute Tm difference (between wild and variant allele)
 * 														when tail is on reverse primer, matching variant allele
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Parameters
 * ──────────────────────────────────────────────────────────────────────────
 * @param {string} targetStrandSeqSnapPrimerRefPoint	Sequence in the snapback primer’s reference orientation (5' -> 3')
 * 														on the strand that shares the same nucleotide pattern as the snapback's
 * 														3' end
 * @param {SNVSite} snvSiteSnapPrimerRefPoint			SNV object in the snapback primer’s reference orientation
 * @param {{start:number,end:number}} bestStemLoc      	Inclusive start/end indices of the best stem (wild type melting temperature
 * 														is closest to the desired melting temperature given by the user) in snapback
 * 														primer reference orientation
 * @param {boolean} tailOnForwardPrimer                 Indicates whether optimal snapback primer has its tail assigned to the forward primer
 *
 * @returns {Promise<MeltingTempDiffs>}                 Object containing absolute melting temperature differences across cases
 *
 * @throws {Error} 										When any of the parameters passed in don't make sense or the API doesn't respond
 * 														correctly
 */
async function calculateMeltingTempDifferences(
	targetStrandSeqSnapPrimerRefPoint,
	snvSiteSnapPrimerRefPoint,
	bestStemLoc,
	tailOnForwardPrimer,
) {
	//──────────────────────────────────────────────────────────────────────────//
	// Parameter Checking									                    //
	//──────────────────────────────────────────────────────────────────────────//

	// 1) Validate sequence (must be non-empty, uppercase A/T/C/G only)
	if (!isValidDNASequence(targetStrandSeqSnapPrimerRefPoint)) {
		throw new Error(
			`Invalid targetStrandSeqSnapPrimerRefPoint: "${targetStrandSeqSnapPrimerRefPoint}". ` +
				'Must be a non-empty, uppercase DNA string containing only A/T/C/G.',
		);
	}

	const SEQ_LEN = targetStrandSeqSnapPrimerRefPoint.length;

	// 2) Validate SNV object with help from helper function
	if (!isValidSNVObject(snvSiteSnapPrimerRefPoint)) {
		throw new Error(
			`Invalid snvSiteSnapPrimerRefPoint: ${JSON.stringify(
				snvSiteSnapPrimerRefPoint,
			)}. Expected { index: number, variantBase: "A"|"T"|"C"|"G" }.`,
		);
	}
	const snvIndex = snvSiteSnapPrimerRefPoint.index;
	if (snvIndex < 0 || snvIndex >= SEQ_LEN) {
		throw new Error(
			`snvSiteSnapPrimerRefPoint.index (${snvIndex}) is out of bounds for ` +
				`sequence length ${SEQ_LEN}.`,
		);
	}

	// 3) Validate bestStemLoc object (shape, only keys, integer bounds, non-empty)
	if (typeof bestStemLoc !== 'object' || bestStemLoc === null) {
		throw new Error('bestStemLoc must be a non-null object.');
	}
	{
		const allowedStemKeys = new Set(['start', 'end']);
		const stemKeys = Object.keys(bestStemLoc);
		for (const k of stemKeys) {
			if (!allowedStemKeys.has(k)) {
				throw new Error(
					`bestStemLoc has unexpected key "${k}". Only {start, end} are allowed.`,
				);
			}
		}
		for (const k of allowedStemKeys) {
			if (!(k in bestStemLoc)) {
				throw new Error(`bestStemLoc missing required key "${k}".`);
			}
		}
	}
	const { start, end } = bestStemLoc;
	if (
		typeof start !== 'number' ||
		!Number.isInteger(start) ||
		!Number.isFinite(start)
	) {
		throw new Error('bestStemLoc.start must be an integer.');
	}
	if (
		typeof end !== 'number' ||
		!Number.isInteger(end) ||
		!Number.isFinite(end)
	) {
		throw new Error('bestStemLoc.end must be an integer.');
	}
	if (start < 0 || end < 0) {
		throw new Error('bestStemLoc.start and bestStemLoc.end must be ≥ 0.');
	}
	if (start > end) {
		throw new Error(
			`bestStemLoc.start (${start}) cannot be greater than bestStemLoc.end (${end}).`,
		);
	}
	if (end >= SEQ_LEN) {
		throw new Error(
			`bestStemLoc.end (${end}) is out of bounds for sequence length ${SEQ_LEN}.`,
		);
	}
	// Non-empty stem (at least one base)
	if (end - start + 1 <= 0) {
		throw new Error('bestStemLoc must span at least one nucleotide.');
	}
	// SNV must lie within the stem for downstream logic/mismatch placement
	if (snvIndex < start || snvIndex > end) {
		throw new Error(
			`snvSiteSnapPrimerRefPoint.index (${snvIndex}) must lie within bestStemLoc ` +
				`[${start}, ${end}] for mismatch positioning.`,
		);
	}

	// 4) Validate tailOnForwardPrimer (boolean)
	if (typeof tailOnForwardPrimer !== 'boolean') {
		throw new Error('tailOnForwardPrimer must be a boolean.');
	}

	//──────────────────────────────────────────────────────────────────────────//
	//								Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1) Creating an object to hold the melting temperature differences
	const meltingTempDiffs = {
		onForwardPrimer: { matchWild: null, matchVariant: null },
		onReversePrimer: { matchWild: null, matchVariant: null },
	};

	// 2) Creating a variable to represent the target strand in the reverse complement
	//	  of the snapback primers reference point
	const targetStrandSeqRevCompSnapPrimerRefPoint = reverseComplement(
		targetStrandSeqSnapPrimerRefPoint,
	);

	// 3) Creating a variable to represent the SNV sithe in the reverse complement
	//	  of the snapback primers reference point
	const snvSiteRevCompSnapPrimerRefPoint = revCompSNV(
		snvSiteSnapPrimerRefPoint,
		targetStrandSeqRevCompSnapPrimerRefPoint.length,
	);

	// 4) Getting the stem location in the reverse complement of the snapback
	//	  primer's reference point.
	const bestStemLocRevCompSnapPrimerRefPoint = {
		start: targetStrandSeqSnapPrimerRefPoint.length - bestStemLoc.end - 1,
		end: targetStrandSeqSnapPrimerRefPoint.length - bestStemLoc.start - 1,
	};

	// 5) Creating the stem sequences any option the stem can be
	const stemSeqWildSamePrimer = targetStrandSeqSnapPrimerRefPoint.slice(
		bestStemLoc.start,
		bestStemLoc.end + 1,
	);
	const stemSeqVariantSamePrimer =
		targetStrandSeqSnapPrimerRefPoint.slice(
			bestStemLoc.start,
			snvSiteSnapPrimerRefPoint.index,
		) +
		snvSiteSnapPrimerRefPoint.variantBase +
		targetStrandSeqSnapPrimerRefPoint.slice(
			snvSiteSnapPrimerRefPoint.index + 1,
			bestStemLoc.end + 1,
		);
	const stemSeqWildRevPrimer = targetStrandSeqRevCompSnapPrimerRefPoint.slice(
		bestStemLocRevCompSnapPrimerRefPoint.start,
		bestStemLocRevCompSnapPrimerRefPoint.end + 1,
	);
	const stemSeqVariantRevPrimer =
		targetStrandSeqRevCompSnapPrimerRefPoint.slice(
			bestStemLocRevCompSnapPrimerRefPoint.start,
			snvSiteRevCompSnapPrimerRefPoint.index,
		) +
		snvSiteRevCompSnapPrimerRefPoint.variantBase +
		targetStrandSeqRevCompSnapPrimerRefPoint.slice(
			snvSiteRevCompSnapPrimerRefPoint.index + 1,
			bestStemLocRevCompSnapPrimerRefPoint.end + 1,
		);

	console.group('Step 5 — Stem sequences at fixed bestStemLoc');
	// Shared context for step 5
	console.log('bestStemLoc, ', bestStemLoc);
	console.log(
		'bestStemLoc (same-primer frame): start=%d, end=%d, len=%d',
		bestStemLoc.start,
		bestStemLoc.end,
		bestStemLoc.end - bestStemLoc.start + 1,
	);
	console.log(
		'bestStemLoc (rev-comp frame): start=%d, end=%d, len=%d',
		bestStemLocRevCompSnapPrimerRefPoint.start,
		bestStemLocRevCompSnapPrimerRefPoint.end,
		bestStemLocRevCompSnapPrimerRefPoint.end -
			bestStemLocRevCompSnapPrimerRefPoint.start +
			1,
	);

	// Same-primer frame (snapback reference orientation)
	console.group('Same primer (snapback reference orientation)');
	console.log('Stem (WILD) sequence: %s', stemSeqWildSamePrimer);
	console.log(
		'  length=%d | global SNV index=%d | SNV base (genomic)=%s | posInStem=%d',
		stemSeqWildSamePrimer.length,
		snvSiteSnapPrimerRefPoint.index,
		targetStrandSeqSnapPrimerRefPoint[snvSiteSnapPrimerRefPoint.index],
		snvSiteSnapPrimerRefPoint.index - bestStemLoc.start,
	);

	console.log('Stem (VARIANT) sequence: %s', stemSeqVariantSamePrimer);
	console.log(
		'  length=%d | global SNV index=%d | SNV base (variant)=%s | posInStem=%d',
		stemSeqVariantSamePrimer.length,
		snvSiteSnapPrimerRefPoint.index,
		snvSiteSnapPrimerRefPoint.variantBase,
		snvSiteSnapPrimerRefPoint.index - bestStemLoc.start,
	);
	console.groupEnd(); // Same primer

	// Reverse-primer frame (reverse complement of snapback reference)
	console.group('Reverse primer (reverse-complement orientation)');
	console.log('Stem (WILD) sequence: %s', stemSeqWildRevPrimer);
	console.log(
		'  length=%d | global SNV index (rev)=%d | SNV base (rev genomic)=%s | posInStem=%d',
		stemSeqWildRevPrimer.length,
		snvSiteRevCompSnapPrimerRefPoint.index,
		targetStrandSeqRevCompSnapPrimerRefPoint[
			snvSiteRevCompSnapPrimerRefPoint.index
		],
		snvSiteRevCompSnapPrimerRefPoint.index -
			bestStemLocRevCompSnapPrimerRefPoint.start,
	);

	console.log('Stem (VARIANT) sequence: %s', stemSeqVariantRevPrimer);
	console.log(
		'  length=%d | global SNV index (rev)=%d | SNV base (rev variant)=%s | posInStem=%d',
		stemSeqVariantRevPrimer.length,
		snvSiteRevCompSnapPrimerRefPoint.index,
		snvSiteRevCompSnapPrimerRefPoint.variantBase,
		snvSiteRevCompSnapPrimerRefPoint.index -
			bestStemLocRevCompSnapPrimerRefPoint.start,
	);
	console.groupEnd(); // Reverse primer

	console.groupEnd(); // Step 5

	// 6) Calculating the loop lengths for each case.
	//    They only depend on which primer the tail is attached to
	const loopLenSamePrimer = getSnapbackLoopLength(
		targetStrandSeqSnapPrimerRefPoint,
		bestStemLoc.start,
	);
	const loopLenRevPrimer = getSnapbackLoopLength(
		targetStrandSeqRevCompSnapPrimerRefPoint,
		bestStemLocRevCompSnapPrimerRefPoint.start,
	);

	// 7) Build mismatch object for wild and variant type on each primer to use in
	//	  snapback Tm calculations

	// Represents the base on the tail matching a wild allele on the
	// REMEMBER: we pass in the stem sequence from the `inner strand`
	//			 This is the strand as seen from `targetStrandSeqSnapPrimerRefPoint`
	//			 or `targetStrandSeqRevCompSnapPrimerRefPoint`.
	//			 The mismatch object should then include the index, relative to the
	//			 start of the stem, and the type which is the base on the tail end
	//			 of the snapback primer, at the location that hybridizes with the SNV.

	const wildTailSamePrimerMismatch = {
		position: snvSiteSnapPrimerRefPoint.index - bestStemLoc.start,
		type: NUCLEOTIDE_COMPLEMENT[
			targetStrandSeqSnapPrimerRefPoint[snvSiteSnapPrimerRefPoint.index]
		],
	};
	const variantTailSamePrimerMismatch = {
		position: snvSiteSnapPrimerRefPoint.index - bestStemLoc.start,
		type: NUCLEOTIDE_COMPLEMENT[snvSiteSnapPrimerRefPoint.variantBase],
	};
	const wildTailRevPrimerMismatch = {
		position:
			snvSiteRevCompSnapPrimerRefPoint.index -
			bestStemLocRevCompSnapPrimerRefPoint.start,
		type: NUCLEOTIDE_COMPLEMENT[
			targetStrandSeqRevCompSnapPrimerRefPoint[
				snvSiteRevCompSnapPrimerRefPoint.index
			]
		],
	};
	const variantTailRevPrimerMismatch = {
		position:
			snvSiteRevCompSnapPrimerRefPoint.index -
			bestStemLocRevCompSnapPrimerRefPoint.start,
		type: NUCLEOTIDE_COMPLEMENT[
			snvSiteRevCompSnapPrimerRefPoint.variantBase
		],
	};

	// 8) Parallelized Tm computations (batch the API calls and then do the math)
	//    - Fires ALL required calculateSnapbackTmWittwer() requests concurrently
	//    - Validates that every call succeeded before computing absolute ΔTm’s

	// 		Kick off ALL eight calls immediately (no awaiting yet):
	const launches = [
		// Same primer, wild stem
		calculateSnapbackTmWittwer(stemSeqWildSamePrimer, loopLenSamePrimer), // 0: same-wild baseline
		calculateSnapbackTmWittwer(
			stemSeqWildSamePrimer,
			loopLenSamePrimer,
			variantTailSamePrimerMismatch,
		), // 1: same-wild + variant tail

		// Same primer, variant stem
		calculateSnapbackTmWittwer(stemSeqVariantSamePrimer, loopLenSamePrimer), // 2: same-variant baseline
		calculateSnapbackTmWittwer(
			stemSeqVariantSamePrimer,
			loopLenSamePrimer,
			wildTailSamePrimerMismatch,
		), // 3: same-variant + wild tail

		// Reverse primer, wild stem
		calculateSnapbackTmWittwer(stemSeqWildRevPrimer, loopLenRevPrimer), // 4: rev-wild baseline
		calculateSnapbackTmWittwer(
			stemSeqWildRevPrimer,
			loopLenRevPrimer,
			variantTailRevPrimerMismatch,
		), // 5: rev-wild + variant tail

		// Reverse primer, variant stem
		calculateSnapbackTmWittwer(stemSeqVariantRevPrimer, loopLenRevPrimer), // 6: rev-variant baseline
		calculateSnapbackTmWittwer(
			stemSeqVariantRevPrimer,
			loopLenRevPrimer,
			wildTailRevPrimerMismatch,
		), // 7: rev-variant + wild tail
	];

	// 		Wait for everything in parallel and verify success:
	const labels = [
		'samePrimer(wild) baseline',
		'samePrimer(wild) + variant-tail mismatch',
		'samePrimer(variant) baseline',
		'samePrimer(variant) + wild-tail mismatch',
		'revPrimer(wild) baseline',
		'revPrimer(wild) + variant-tail mismatch',
		'revPrimer(variant) baseline',
		'revPrimer(variant) + wild-tail mismatch',
	];

	const settled = await Promise.allSettled(launches);
	const failIdx = settled.findIndex((r) => r.status === 'rejected');
	if (failIdx !== -1) {
		const reason = settled[failIdx].reason;
		// Surface a precise, actionable error (bubbles to your existing try/catch)
		throw new Error(
			`calculateSnapbackTmWittwer failed for ${labels[failIdx]}: ${
				reason?.message ?? String(reason)
			}`,
		);
	}

	// 		Unpack numeric results in the same order as launched:
	const [
		sameWild_base,
		sameWild_withVariantTail,
		sameVar_base,
		sameVar_withWildTail,
		revWild_base,
		revWild_withVariantTail,
		revVar_base,
		revVar_withWildTail,
	] = settled.map((r) => /** @type {number} */ (r.value));

	// 9) Compute absolute Tm differences with API results
	const meltingTempDiffSamePrimerMatchWild = Math.abs(
		sameWild_base - sameVar_withWildTail,
	); // remember the tail is what we choose and what stays the same. What it anneals to causese the melting temperatature differences.
	const meltingTempDiffSamePrimerMatchVariant = Math.abs(
		sameVar_base - sameWild_withVariantTail,
	);
	const meltingTempDiffRevPrimerMatchWild = Math.abs(
		revWild_base - revVar_withWildTail,
	);
	const meltingTempDiffRevPrimerMatchVariant = Math.abs(
		revVar_base - revWild_withVariantTail,
	);

	// 10) Map Tm difference values to their correct orientation depending on `tailOnForwardPrimer`.)
	if (tailOnForwardPrimer) {
		meltingTempDiffs.onForwardPrimer.matchWild =
			meltingTempDiffSamePrimerMatchWild;
		meltingTempDiffs.onForwardPrimer.matchVariant =
			meltingTempDiffSamePrimerMatchVariant;
		meltingTempDiffs.onReversePrimer.matchWild =
			meltingTempDiffRevPrimerMatchWild;
		meltingTempDiffs.onReversePrimer.matchVariant =
			meltingTempDiffRevPrimerMatchVariant;
	} else {
		meltingTempDiffs.onForwardPrimer.matchWild =
			meltingTempDiffRevPrimerMatchWild;
		meltingTempDiffs.onForwardPrimer.matchVariant =
			meltingTempDiffRevPrimerMatchVariant;
		meltingTempDiffs.onReversePrimer.matchWild =
			meltingTempDiffSamePrimerMatchWild;
		meltingTempDiffs.onReversePrimer.matchVariant =
			meltingTempDiffSamePrimerMatchVariant;
	}

	// 11) Return the results
	return meltingTempDiffs;
}

/**
 * Decide whether the snapback tail should be appended to the forward primer
 * (target-sequence strand) or to the reverse primer, and which base at the SNV
 * position, on the snapback tail, maximizes the absolute melting-temperature
 * difference between the wild-type and variant stems in the initial stem region
 * (length 2 × SNV_BASE_BUFFER + 1).
 *
 * Assumptions:
 * - targetSeqStrand is a valid uppercase DNA string (A, T, C, G) written 5'→3'.
 * - snvSite passes isValidSNVObject and its index is at least SNV_BASE_BUFFER
 *   bases from both sequence ends.
 *
 * Parameters:
 * @param {string}   targetSeqStrand  The strand to which the forward primer binds.
 * @param {SNVSite}  snvSite          { index: number, variantBase: "A"|"T"|"C"|"G" }
 *
 * @returns {Promise<{
 *   tailOnForwardPrimer   : boolean,
 *   snapbackTailBaseAtSNV     : string,
 *   snapbackTailMatchesWild: boolean
 * }>}
 *
 * @throws {Error} If inputs are malformed or violate positional constraints.
 */
async function useForwardPrimer(targetSeqStrand, snvSite) {
	//──────────────────────────────────────────────────────────────────────────//
	// Parameter checking                                                      //
	//──────────────────────────────────────────────────────────────────────────//

	// 1. Validate targetSeqStrand
	if (!isValidDNASequence(targetSeqStrand)) {
		throw new Error(
			'targetSeqStrand must be a non-empty uppercase DNA string containing only A, T, C, or G.',
		);
	}

	// 2. Validate snvSite structure and content
	if (!isValidSNVObject(snvSite)) {
		throw new Error(
			`snvSite is invalid: ${JSON.stringify(
				snvSite,
			)}. Expected { index: non-negative integer, variantBase: "A"|"T"|"C"|"G" }.`,
		);
	}

	// 3. Ensure snvSite.index is within sequence bounds
	if (snvSite.index >= targetSeqStrand.length) {
		throw new Error(
			`snvSite.index (${snvSite.index}) exceeds sequence length ${targetSeqStrand.length}.`,
		);
	}

	// 4. Ensure the SNV is sufficiently distant from both sequence ends
	if (
		snvSite.index < SNV_BASE_BUFFER ||
		snvSite.index > targetSeqStrand.length - SNV_BASE_BUFFER - 1
	) {
		throw new Error(
			`SNV at index ${snvSite.index} is too close to a sequence end; ` +
				`need at least ${SNV_BASE_BUFFER} perfectly matched bases flanking it.`,
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

	// 4) Slice out the "init stem" region from each strand. Target corresponds to the
	//	  target
	const targetInitStem = targetSeqStrand.slice(
		initStemLoc.start,
		initStemLoc.end + 1,
	);
	const compInitStem = revCompTargetSeqStrand.slice(
		compInitStemLoc.start,
		compInitStemLoc.end + 1,
	);

	// 6) SNV is at position {SNV_BASE_BUFFER} in these 2*{SNV_BASE_BUFFER}+1 slices
	const mismatchPos = SNV_BASE_BUFFER;

	// 7) Evaluate Tm differences for snapback tail on target strand
	const tailOnForwardPrimerScenario =
		await evaluateSnapbackTailMatchingOptions(
			targetInitStem,
			mismatchPos,
			snvSite.variantBase,
		);

	// 8) Evaluate Tm differences for snapback tail on complementary strand
	const tailOnReversePrimerScenario =
		await evaluateSnapbackTailMatchingOptions(
			compInitStem,
			mismatchPos,
			revCompSnvSite.variantBase,
		);

	console.log('tailOnForwardPrimerScenario', tailOnForwardPrimerScenario);
	console.log('tailOnReversePrimerScenario', tailOnReversePrimerScenario);

	// 0) Compare which scenario yields the bigger Tm difference
	if (
		tailOnForwardPrimerScenario.bestTmDifference >
		tailOnReversePrimerScenario.bestTmDifference
	) {
		return {
			tailOnForwardPrimer: true,
			bestSnapbackTailBaseAtSNV:
				tailOnForwardPrimerScenario.bestSnapbackTailBaseAtSNV,
			snapbackTailMatchesWild:
				tailOnForwardPrimerScenario.snapbackTailMatchesWild,
		};
	} else {
		return {
			tailOnForwardPrimer: false,
			bestSnapbackTailBaseAtSNV:
				tailOnReversePrimerScenario.bestSnapbackTailBaseAtSNV,
			snapbackTailMatchesWild:
				tailOnReversePrimerScenario.snapbackTailMatchesWild,
		};
	}
}

/**
 * Evaluates which snapback-tail base (wild-matching vs. variant-matching)
 * maximises the Tm difference between wild-type and variant stems in the
 * initial “seed” slice.
 *
 * Two cases are compared:
 *  1. Tail matches the wild base → variant stem contains the mismatch.
 *  2. Tail matches the variant base → wild stem contains the mismatch.
 *
 * The scenario with the larger |ΔTm| wins.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Parameters, Returns, and Errors
 * ──────────────────────────────────────────────────────────────────────────
 * @param {string} initStem      Slice of length (2·SNV_BASE_BUFFER + 1) with
 *                               the wild base at `mismatchPos`.
 * @param {number} mismatchPos   Index of the SNV within `initStem`
 *                               (normally SNV_BASE_BUFFER).
 * @param {string} variantBase   Variant base as it appears on the target
 *                               strand (A/T/C/G).
 *
 * @returns {Promise<{
 *   bestSnapbackTailBaseAtSNV        	: string,  // Base to place in the tail
 *   bestDifference          			: number,  // Larger |ΔTm| in °C
 *   snapbackTailMatchesWild 			: boolean  // true → tail matches wild base
 * }>}
 *
 * @throws {Error} If any argument is invalid.
 */
async function evaluateSnapbackTailMatchingOptions(
	initStem,
	mismatchPos,
	variantBase,
) {
	//──────────────────────────────────────────────────────────────────────────//
	//							Parameter Checking								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1. initStem must be a valid DNA sequence
	if (!isValidDNASequence(initStem)) {
		throw new Error(`initStem must be a non-empty A/T/C/G string.`);
	}

	// 2. mismatchPos must be a valid index inside initStem
	if (
		typeof mismatchPos !== 'number' ||
		!Number.isInteger(mismatchPos) ||
		mismatchPos < 0 ||
		mismatchPos >= initStem.length
	) {
		throw new Error(
			`mismatchPos (${mismatchPos}) must be an integer between 0 and ${
				initStem.length - 1
			}.`,
		);
	}

	// 3. variantBase must be a single valid nucleotide
	if (
		typeof variantBase !== 'string' ||
		variantBase.length !== 1 ||
		!VALID_BASES.has(variantBase)
	) {
		throw new Error(
			`variantBase must be one character A, T, C, or G. Received "${variantBase}".`,
		);
	}

	// 4. variantBase must differ from the wild base at mismatchPos
	const wildBase = initStem[mismatchPos];
	if (variantBase === wildBase) {
		throw new Error(
			`variantBase must differ from wild base "${wildBase}" at mismatchPos.`,
		);
	}

	//──────────────────────────────────────────────────────────────────────────//
	//								Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1) Get Tm for a stem with the wild base where the snapback tail matches it
	const wildMatchTmPromise = await getOligoTm(initStem);

	// 2) Get Tm for a stem with the variant allele where the snapback tail matches it
	const variantInitStem =
		initStem.slice(0, mismatchPos) +
		variantBase +
		initStem.slice(mismatchPos + 1);

	const variantMatchTmPromise = await getOligoTm(variantInitStem);

	// 3) Resolve both promises
	//    I did this to learn more about aynchronous code, I could do it more
	//    throughout my project, but it won't speed things up much.
	const [wildMatchTm, variantMatchTm] = await Promise.all([
		wildMatchTmPromise,
		variantMatchTmPromise,
	]);

	// 2) Scenario A: Wild-matching snapback tail is mismatched when forward primer
	//    part anneals to complement sequence with variant type base
	const wildMatchingSnapbackTailToVariantMismatchObj = {
		position: mismatchPos,
		// At the mismatch location, the snapback tail matches the wild type base
		type: NUCLEOTIDE_COMPLEMENT[wildBase],
	};
	const wildMatchingSnapbackTailToVariantTm = await getOligoTm(
		variantInitStem,
		wildMatchingSnapbackTailToVariantMismatchObj,
	);
	const wildMatchingSnapbackTailTmDiff = Math.abs(
		wildMatchTm - wildMatchingSnapbackTailToVariantTm,
	);

	// 3) Scenario B: Variant-matching snapback tail
	//    (mis)matches when forward primer anneals to complement sequend with wild type base
	const variantMatchingSnapbackTailToWildMismatchObj = {
		position: mismatchPos,
		// At the mismatch, the snapback nucleotide will be the complement of the variant type base
		type: NUCLEOTIDE_COMPLEMENT[variantBase],
	};
	const variantMatchingSnapbackTailToWildTm = await getOligoTm(
		initStem,
		variantMatchingSnapbackTailToWildMismatchObj,
	);
	const variantMatchingSnapbackTailTmDiff = Math.abs(
		variantMatchTm - variantMatchingSnapbackTailToWildTm,
	);

	/* ───────────── DEBUG LOG BLOCK (paste before the return) ───────────── */

	console.group('evaluateSnapbackTailMatchingOptions — debug');

	// Inputs & slices
	console.log('initStem:', initStem, 'len=', initStem.length);
	console.log('mismatchPos:', mismatchPos);
	console.log('wildBase @mismatchPos:', initStem[mismatchPos]);
	console.log('variantBase:', variantBase);
	console.log('variantInitStem:', variantInitStem);

	// Base complements (for sanity)
	console.log('comp(wildBase):', NUCLEOTIDE_COMPLEMENT[wildBase]);
	console.log('comp(variantBase):', NUCLEOTIDE_COMPLEMENT[variantBase]);

	// Top-strand Tms (no mismatch)
	console.log('wildMatchTm (initStem):', wildMatchTm);
	console.log('variantMatchTm (variantInitStem):', variantMatchTm);

	// Scenario A (tail matches WILD → mismatch vs VARIANT top)
	console.log(
		'Scenario A mismatch obj:',
		wildMatchingSnapbackTailToVariantMismatchObj,
	);
	console.log(
		'Scenario A Tm (variantInitStem + A_mismatch):',
		wildMatchingSnapbackTailToVariantTm,
	);
	console.log(
		'Scenario A ΔTm = |wildMatchTm - A_Tm|:',
		Math.abs(wildMatchTm - wildMatchingSnapbackTailToVariantTm),
		'→ stored as wildMatchingSnapbackTailTmDiff=',
		wildMatchingSnapbackTailTmDiff,
	);

	// Scenario B (tail matches VARIANT → mismatch vs WILD top)
	console.log(
		'Scenario B mismatch obj:',
		variantMatchingSnapbackTailToWildMismatchObj,
	);
	console.log(
		'Scenario B Tm (initStem + B_mismatch):',
		variantMatchingSnapbackTailToWildTm,
	);
	console.log(
		'Scenario B ΔTm = |variantMatchTm - B_Tm|:',
		Math.abs(variantMatchTm - variantMatchingSnapbackTailToWildTm),
		'→ stored as variantMatchingSnapbackTailTmDiff=',
		variantMatchingSnapbackTailTmDiff,
	);

	// Summary table
	console.table([
		{
			Scenario: 'A: tail matches WILD',
			topStrand: 'variantInitStem',
			mismatchType: wildMatchingSnapbackTailToVariantMismatchObj.type,
			BasePlacedInTail: NUCLEOTIDE_COMPLEMENT[wildBase],
			NoMismatchTm: wildMatchTm,
			ScenarioTm: wildMatchingSnapbackTailToVariantTm,
			DeltaTm: wildMatchingSnapbackTailTmDiff,
		},
		{
			Scenario: 'B: tail matches VARIANT',
			topStrand: 'initStem',
			mismatchType: variantMatchingSnapbackTailToWildMismatchObj.type,
			BasePlacedInTail: NUCLEOTIDE_COMPLEMENT[variantBase],
			NoMismatchTm: variantMatchTm,
			ScenarioTm: variantMatchingSnapbackTailToWildTm,
			DeltaTm: variantMatchingSnapbackTailTmDiff,
		},
	]);

	// Decision preview (before actual return)
	const _chooseA =
		wildMatchingSnapbackTailTmDiff > variantMatchingSnapbackTailTmDiff;
	console.log('Decision preview:', {
		chooseA: _chooseA,
		reason: _chooseA
			? 'Scenario A ΔTm > Scenario B ΔTm'
			: 'Scenario B ΔTm ≥ Scenario A ΔTm',
		scenarioA_DeltaTm: wildMatchingSnapbackTailTmDiff,
		scenarioB_DeltaTm: variantMatchingSnapbackTailTmDiff,
		chosenTailBase: _chooseA
			? NUCLEOTIDE_COMPLEMENT[wildBase]
			: NUCLEOTIDE_COMPLEMENT[variantBase],
		snapbackTailMatchesWildPreview: _chooseA,
	});

	// Tie/suspicious conditions
	if (wildMatchingSnapbackTailTmDiff === variantMatchingSnapbackTailTmDiff) {
		console.warn('ΔTm tie detected: A == B. Check inputs/orientation.');
	}
	if (variantBase === wildBase) {
		console.warn(
			'variantBase equals wildBase at mismatchPos (unexpected).',
		);
	}

	console.groupEnd();
	/* ───────────── END DEBUG LOG BLOCK ───────────── */

	// 4) Pick whichever scenario yields the larger difference
	if (wildMatchingSnapbackTailTmDiff > variantMatchingSnapbackTailTmDiff) {
		// Scenario A wins
		// Snapback tail should match wild
		return {
			bestSnapbackTailBaseAtSNV: NUCLEOTIDE_COMPLEMENT[wildBase],
			bestTmDifference: wildMatchingSnapbackTailTmDiff,
			snapbackTailMatchesWild: true,
		};
	} else {
		// Scenario B wins
		// Snapback tail should match variant
		return {
			bestSnapbackTailBaseAtSNV: NUCLEOTIDE_COMPLEMENT[variantBase],
			bestTmDifference: variantMatchingSnapbackTailTmDiff,
			snapbackTailMatchesWild: false,
		};
	}
}

/**
 * Retrieves the melting temperature (Tm) of a perfectly matched or
 * single-mismatch DNA stem by querying the dna-utah.org Santa Lucia CGI.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Assumptions
 * ──────────────────────────────────────────────────────────────────────────
 * – `seq` is an uppercase DNA string (A/T/C/G) validated by `isValidDNASequence`.
 * – If `mismatch` is supplied it must pass `isValidMismatchObject`.
 * – Global constants (MG, MONO, etc.) are defined in module scope.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Type definitions
 * ──────────────────────────────────────────────────────────────────────────
 * @typedef {Object} Mismatch
 * @property {number} position   Zero-based index within `seq`
 * @property {string} type       Intended base on the opposite strand (A/T/C/G)
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Parameters, Returns, and Errors
 * ──────────────────────────────────────────────────────────────────────────
 * @param  {string}    seq              Fully matched reference sequence (5'→3')
 * @param  {Mismatch} [mismatch]        Optional mismatch specification
 *
 * @returns {Promise<number>}           Melting temperature (°C) rounded to
 *                                      `TM_DECIMAL_PLACES`
 *
 * @throws {Error}                      If inputs are invalid, the network
 *                                      request fails, or the CGI response
 *                                      cannot be parsed
 */
async function getOligoTm(seq, mismatch) {
	//──────────────────────────────────────────────────────────────────────//
	// Parameter Checking                                                   //
	//──────────────────────────────────────────────────────────────────────//

	// 1. Validate the DNA sequence
	if (!isValidDNASequence(seq)) {
		throw new Error(
			`Invalid DNA sequence: "${seq}". Must be non-empty and contain only A, T, C, or G.`,
		);
	}

	// 2. Validate the mismatch object (if provided)
	if (mismatch !== undefined && mismatch !== null) {
		// 2.a  Shape and content
		if (!isValidMismatchObject(mismatch)) {
			throw new Error(
				`Invalid mismatch object: ${JSON.stringify(
					mismatch,
				)}. Expected { position: int, type: "A"|"T"|"C"|"G" }.`,
			);
		}
		// 2.b  Position must lie within sequence bounds
		if (mismatch.position >= seq.length) {
			throw new Error(
				`Mismatch position (${mismatch.position}) exceeds sequence length ${seq.length}.`,
			);
		}
	}

	//──────────────────────────────────────────────────────────────────────────//
	// Function Logic                                                          //
	//──────────────────────────────────────────────────────────────────────────//

	// 1. Build the mmseq string if a mismatch is present
	let mmSeq = null;
	if (mismatch) {
		mmSeq = buildMismatchSequenceForAPI(seq, mismatch);
	}

	// 2. Assemble the query URL
	let apiURL = API_URL;
	apiURL += `?mg=${MG}`;
	apiURL += `&mono=${MONO}`;
	apiURL += `&seq=${seq.toLowerCase()}`;
	apiURL += `&tparam=${T_PARAM}`;
	apiURL += `&saltcalctype=${SALT_CALC_TYPE}`;
	apiURL += `&otype=${O_TYPE}`;

	//apiURL += `&concentration=${CONC}`;
	//apiURL += `&limitingconc=${LIMITING_CONC}`;
	apiURL += `&decimalplaces=${TM_DECIMAL_PLACES}`;
	if (USE_TOKEN && API_TOKEN) {
		apiURL += `&token=${API_TOKEN}`;
	}
	if (mmSeq) {
		apiURL += `&mmseq=${mmSeq}`;
	}

	// 3. If proxying, encode the target and prepend the proxy URL
	const finalURL = USE_PROXY
		? `${PROXY_URL}?url=${encodeURIComponent(apiURL)}`
		: apiURL;

	// 4. Fetch the response
	const res = await fetch(finalURL);
	if (!res.ok) {
		throw new Error(`Network error: ${res.status} – ${res.statusText}`);
	}
	const rawHtml = await res.text();

	// 5. Extract the Tm (wild-type <tm> or mismatch <mmtm>)
	const tmVal = parseTmFromResponse(rawHtml, Boolean(mismatch));

	// 6. Validate that a numeric Tm was found
	if (tmVal === null) {
		throw new Error('Tm value not found or unparsable in server response.');
	}

	// 7. Return the temperature
	return tmVal;
}

/**
 * Retrieves the melting temperature (Tm) of a primer by querying the
 * dna-utah.org Santa Lucia CGI with otype=primer.
 *
 * This mirrors getOligoTm but targets the primer-specific calculation mode.
 *
 * @param  {string}    seq              Primer sequence (5'→3')
 * @returns {Promise<number>}           Melting temperature (°C) rounded to
 *                                      `TM_DECIMAL_PLACES`
 * @throws {Error}                      If inputs are invalid, the network
 *                                      request fails, or the CGI response
 *                                      cannot be parsed
 */
async function getPrimerTm(seq) {
	//──────────────────────────────────────────────────────────────────────//
	// Parameter Checking                                                   //
	//──────────────────────────────────────────────────────────────────────//

	// 1. Validate the DNA sequence
	if (!isValidDNASequence(seq)) {
		throw new Error(
			`Invalid DNA sequence: "${seq}". Must be non-empty and contain only A, T, C, or G.`,
		);
	}

	//──────────────────────────────────────────────────────────────────────────//
	// Function Logic                                                          //
	//──────────────────────────────────────────────────────────────────────────//

	// 1. Assemble the query URL
	let apiURL = API_URL;
	apiURL += `?mg=${MG}`; // 2.2 currently
	apiURL += `&mono=${MONO}`; // 20 currently
	apiURL += `&seq=${seq.toLowerCase()}`;
	apiURL += `&tparam=${T_PARAM}`; //SantaLuciaHicks currently
	apiURL += `&saltcalctype=${SALT_CALC_TYPE}`; // 'bpdenominator' currenltly
	apiURL += `&otype=${PRIMER_O_TYPE}`; // 'primer' currently

	//apiURL += `&concentration=${CONC}`;
	//apiURL += `&limitingconc=${LIMITING_CONC}`;
	apiURL += `&decimalplaces=${TM_DECIMAL_PLACES}`;
	if (USE_TOKEN && API_TOKEN) {
		apiURL += `&token=${API_TOKEN}`;
	}

	// 2. If proxying, encode the target and prepend the proxy URL
	const finalURL = USE_PROXY
		? `${PROXY_URL}?url=${encodeURIComponent(apiURL)}`
		: apiURL;

	// 3. Fetch the response
	const res = await fetch(finalURL);
	if (!res.ok) {
		throw new Error(`Network error: ${res.status} – ${res.statusText}`);
	}
	const rawHtml = await res.text();

	// 4. Extract the Tm
	const tmVal = parseTmFromResponse(rawHtml);

	// 5. Validate that a numeric Tm was found
	if (tmVal === null) {
		throw new Error('Tm value not found or unparsable in server response.');
	}

	// 6. Return the temperature
	return tmVal;
}

/**
 * Retrieves ΔH°, ΔS°, and the salt-correction term for a DNA duplex by querying
 * the dna-utah Tm API. Inputs for concentrations are in µM.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Assumptions
 * ──────────────────────────────────────────────────────────────────────────
 * - `seq` is a valid uppercase DNA string (A/T/C/G) validated via isValidDNASequence.
 * - If `mismatch` is supplied it must pass isValidMismatchObject and be in-bounds.
 * - Global constants (MG, MONO, T_PARAM, SALT_CALC_TYPE, O_TYPE, TM_DECIMAL_PLACES,
 *   API_URL, API_TOKEN, USE_PROXY, PROXY_URL) are defined in module scope.
 * - The endpoint returns:
 *     <dH> cal/mol </dH>
 *     <dS> cal/K/mol </dS>
 *     <saltCorrection> °C </saltCorrection>
 *   We convert dH → kcal/mol for convenience/consistency with your Tm code.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Parameters, Returns, and Errors
 * ──────────────────────────────────────────────────────────────────────────
 * @typedef {Object} Mismatch
 * @property {number} 	position   			Zero-based index within `seq`
 * @property {string} 	type       			Intended base on the opposite strand (A/T/C/G)
 *
 * @param  {string}   	seq               	DNA sequence 5'→3' (uppercase A/T/C/G)
 * @param  {number}   	concentration     	Main strand concentration in µM (> 0)
 * @param  {number}   	limitingConc      	Limiting strand concentration in µM (> 0)
 * @param  {Mismatch} 	[mismatch]        	Optional mismatch specification
 *
 * @returns {Promise<{ dH:number, dS:number, saltCorrection:number }>}
 *          dH in kcal/mol, dS in cal/K/mol, saltCorrection in °C.
 *
 * @throws  {Error} 						on invalid params, network failure,
 * 											or unparsable response.
 */
async function getThermoParams(seq, concentration, limitingConc, mismatch) {
	//──────────────────────────────────────────────────────────────────────//
	// Parameter Checking                                                   //
	//──────────────────────────────────────────────────────────────────────//

	// 1) Sequence
	if (!isValidDNASequence(seq)) {
		throw new Error(
			`Invalid DNA sequence: "${seq}". Must be non-empty and contain only A, T, C, or G.`,
		);
	}

	// 2) Concentrations (µM) — optional. Validate only if provided.
	if (concentration != null) {
		if (
			typeof concentration !== 'number' ||
			!Number.isFinite(concentration) ||
			concentration <= 0
		) {
			throw new Error(
				`concentration must be a positive, finite number in µM when provided. Received: ${concentration}`,
			);
		}
	}
	if (limitingConc != null) {
		if (
			typeof limitingConc !== 'number' ||
			!Number.isFinite(limitingConc) ||
			limitingConc <= 0
		) {
			throw new Error(
				`limitingConc must be a positive, finite number in µM when provided. Received: ${limitingConc}`,
			);
		}
	}

	// 3) Mismatch (optional)
	if (mismatch !== undefined && mismatch !== null) {
		if (!isValidMismatchObject(mismatch)) {
			throw new Error(
				`Invalid mismatch object: ${JSON.stringify(
					mismatch,
				)}. Expected { position: int, type: "A"|"T"|"C"|"G" }.`,
			);
		}
		if (mismatch.position >= seq.length) {
			throw new Error(
				`Mismatch position (${mismatch.position}) exceeds sequence length ${seq.length}.`,
			);
		}
	}

	//──────────────────────────────────────────────────────────────────────//
	// Request construction                                                //
	//──────────────────────────────────────────────────────────────────────//

	// 1. Build the mmseq string if a mismatch is present
	let mmSeq = null;
	if (mismatch) {
		mmSeq = buildMismatchSequenceForAPI(seq, mismatch);
	}

	// 2. Assemble the query URL
	let apiURL = API_URL;
	apiURL += `?mg=${MG}`;
	apiURL += `&mono=${MONO}`;
	apiURL += `&seq=${seq.toLowerCase()}`;
	apiURL += `&tparam=${T_PARAM}`;
	apiURL += `&saltcalctype=${SALT_CALC_TYPE}`;
	apiURL += `&concentration=${CONC}`;
	apiURL += `&limitingconc=${LIMITING_CONC}`;
	apiURL += `&otype=${O_TYPE}`;
	apiURL += `&decimalplaces=${TM_DECIMAL_PLACES}`;
	if (USE_TOKEN && API_TOKEN) {
		apiURL += `&token=${API_TOKEN}`;
	}
	if (mmSeq) {
		apiURL += `&mmseq=${mmSeq}`;
	}

	// 3. If proxying, encode the target and prepend the proxy URL
	const finalURL = USE_PROXY
		? `${PROXY_URL}?url=${encodeURIComponent(apiURL)}`
		: apiURL;

	// 4. Fetch the response
	const res = await fetch(finalURL);
	if (!res.ok) {
		throw new Error(`Network error: ${res.status} – ${res.statusText}`);
	}
	const rawHtml = await res.text();

	// 5. Parse thermo parameters and return them  via the dedicated parser
	const parsedThermoParams = parseThermoParamsFromResponse(rawHtml);
	return parsedThermoParams;
}

/**
 * Constructs a snapback stem on the selected primer-bearing strand so that the
 * wild-type melting temperature (Tm) approaches `targetSnapMeltTemp`.  Growth
 * proceeds symmetrically (right, then left, repeating) while:
 *   • Maintaining ≥ SNV_BASE_BUFFER perfectly matched bases between the SNV
 *     and each primer binding site.
 *   • Keeping the loop as short as possible, with one extra loop-side
 *     mismatch only when the primer 5' base complements the base immediately
 *     left of the stem.
 *   • Leaving the SNV centred (or as near-centred as sequence boundaries allow)
 *     within the final stem.
 *
 * Extension stops when either primer boundaries are reached or the computed
 * wild-type Tm meets/exceeds `targetSnapMeltTemp`.  If the resulting Tm is
 * still below MINIMUM_TARGET_SNAPBACK_MELTING_TEMP an error is thrown.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Assumptions
 * ──────────────────────────────────────────────────────────────────────────
 * – `targetStrandSeqSnapPrimerRefPoint` is a valid uppercase DNA string (5'→3').
 * – `primerLensSnapPrimerRefPoint.primerLen` and `.compPrimerLen` ≥ MIN_PRIMER_LEN.
 * – The SNV lies ≥ SNV_BASE_BUFFER bases away from both primers.
 * – `snapbackTailBaseAtSNV` is the complement of either the wild or variant base
 *   (chosen earlier for maximal |ΔTm| in the seed stem).
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Type definitions
 * ──────────────────────────────────────────────────────────────────────────
 * @typedef {Object} SNVSiteRefPoint
 *         @property {number} index        0-based SNV position on this strand
 *         @property {string} variantBase  "A" | "T" | "C" | "G"
 *
 * @typedef {Object} PrimerLensRefPoint
 *         @property {number} primerLen      Length of primer on this strand
 *         @property {number} compPrimerLen  Length of complementary primer
 *
 * @typedef {Object} MeltingTemp
 *         @property {number} wildTm     Snapback Tm on wild-type allele (°C)
 *         @property {number} variantTm  Snapback Tm on variant allele (°C)
 *
 * @typedef {Object} StemLoc
 *         @property {number} start  Inclusive 0-based start index of stem
 *         @property {number} end    Inclusive 0-based end   index of stem
 *
 * @typedef {Object} CreateStemReturn
 *         @property {StemLoc}   bestStemLoc               Finalised stem location
 *         @property {MeltingTemp} meltingTemps        Wild / variant Tm values
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Parameters, Returns, and Errors
 * ──────────────────────────────────────────────────────────────────────────
 * @param {string}              targetStrandSeqSnapPrimerRefPoint
 *                                  DNA sequence (5'→3') of the strand that
 *                                  will receive the snapback tail.
 * @param {SNVSiteRefPoint}     snvSiteSnapPrimerRefPoint
 *                                  SNV description in this strand’s coordinates.
 * @param {PrimerLensRefPoint}  primerLensSnapPrimerRefPoint
 *                                  Object holding `primerLen` and `compPrimerLen`.
 * @param {string}              snapbackTailBaseAtSNV
 *                                  Complement base inserted at SNV in the tail.
 * @param {boolean}             matchesWild
 *                                  true → tail complements wild allele,
 *                                  false → tail complements variant allele.
 * @param {number}              targetSnapMeltTemp
 *                                  Desired wild-type stem Tm (°C).
 *
 * @returns {CreateStemReturn}  Object containing stem location and the Tm data.
 *
 * @throws {Error}  If any argument is malformed, the SNV is too close to a
 *                  primer, or a stem meeting temperature/length constraints
 *                  cannot be constructed.
 */
async function createStem(
	targetStrandSeqSnapPrimerRefPoint,
	snvSiteSnapPrimerRefPoint,
	primerLensSnapPrimerRefPoint,
	snapbackTailBaseAtSNV,
	matchesWild,
	targetSnapMeltTemp,
) {
	//──────────────────────────────────────────────────────────────────────────//
	//                          Parameter Checking                              //
	//──────────────────────────────────────────────────────────────────────────//

	// 1. Validate the target strand sequence
	if (!isValidDNASequence(targetStrandSeqSnapPrimerRefPoint)) {
		throw new Error(
			'Invalid targetStrandSeqSnapPrimerRefPoint: must be a non-empty A/T/C/G string.',
		);
	}

	// 2. Validate the SNV object
	if (!isValidSNVObject(snvSiteSnapPrimerRefPoint)) {
		throw new Error(
			`snvSiteSnapPrimerRefPoint is invalid: ${JSON.stringify(
				snvSiteSnapPrimerRefPoint,
			)}`,
		);
	}
	// 2a. Ensure SNV index is within sequence bounds
	if (
		snvSiteSnapPrimerRefPoint.index >=
		targetStrandSeqSnapPrimerRefPoint.length
	) {
		throw new Error(
			`snvSiteSnapPrimerRefPoint.index (${snvSiteSnapPrimerRefPoint.index}) exceeds sequence length ${targetStrandSeqSnapPrimerRefPoint.length}.`,
		);
	}

	// 3. Validate snapbackTailBaseAtSNV
	if (
		typeof snapbackTailBaseAtSNV !== 'string' ||
		snapbackTailBaseAtSNV.length !== 1 ||
		!VALID_BASES.has(snapbackTailBaseAtSNV)
	) {
		throw new Error(
			`snapbackTailBaseAtSNV ("${snapbackTailBaseAtSNV}") must be one of "A", "T", "C", or "G".`,
		);
	}

	// 4. Validate matchesWild flag
	if (typeof matchesWild !== 'boolean') {
		throw new Error('matchesWild must be a boolean.');
	}

	// 5. Validate primerLensSnapPrimerRefPoint structure
	if (
		typeof primerLensSnapPrimerRefPoint !== 'object' ||
		primerLensSnapPrimerRefPoint === null ||
		Array.isArray(primerLensSnapPrimerRefPoint) ||
		!('primerLen' in primerLensSnapPrimerRefPoint) ||
		!('compPrimerLen' in primerLensSnapPrimerRefPoint)
	) {
		throw new Error(
			'primerLensSnapPrimerRefPoint must be an object with integer properties "primerLen" and "compPrimerLen".',
		);
	}

	const { primerLen, compPrimerLen } = primerLensSnapPrimerRefPoint;
	for (const [name, len] of [
		['primerLen', primerLen],
		['compPrimerLen', compPrimerLen],
	]) {
		if (
			typeof len !== 'number' ||
			!Number.isInteger(len) ||
			len < MIN_PRIMER_LEN
		) {
			throw new Error(
				`${name} must be an integer ≥ ${MIN_PRIMER_LEN}. Received ${len}.`,
			);
		}
	}

	// 6. Ensure primer regions fit within the sequence
	const seqLen = targetStrandSeqSnapPrimerRefPoint.length;
	if (primerLen + compPrimerLen >= seqLen) {
		throw new Error(
			`primerLen (${primerLen}) + compPrimerLen (${compPrimerLen}) cannot equal or exceed sequence length (${seqLen}).`,
		);
	}

	// 7. Ensure the SNV is sufficiently distant from both primers
	if (
		snvTooCloseToPrimer(
			snvSiteSnapPrimerRefPoint.index,
			primerLen,
			compPrimerLen,
			seqLen,
		)
	) {
		throw new Error(
			`SNV at index ${snvSiteSnapPrimerRefPoint.index} is within ${SNV_BASE_BUFFER} bases of a primer binding site.`,
		);
	}

	// 8. Validate the targetSnapMeltTemp
	if (
		typeof targetSnapMeltTemp !== 'number' ||
		!Number.isFinite(targetSnapMeltTemp) ||
		targetSnapMeltTemp <= 0
	) {
		throw new Error(
			`targetSnapMeltTemp must be a positive, finite number. Received ${targetSnapMeltTemp}.`,
		);
	}

	//──────────────────────────────────────────────────────────────────────────//
	//								Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//

	//// We will keep enlarging the stem, right then left... (as long as we are not up against the primers), keeping
	//// track of the melting temperature of the snapback for the wild-type allele, until we go over the desired meling
	//// temperature or we run out of viable stem location

	// 1. Initialize variable to hold the snapback melting temperature for the wild type allele that is closest to the desired
	// snapback melting temperature for the wild type allele. Also initialize a variable for the corresponding melting
	// temperature of the variant snapback melting temperature. Finally, initialize the variable for the corresponding stem locations
	let bestWildTm = null;
	let correspondingVariantStemTm = null;
	let bestStemLoc = { start: null, end: null };

	// 2. Initialize stem region
	const snvIndex = snvSiteSnapPrimerRefPoint.index;
	let stemStart = snvIndex - SNV_BASE_BUFFER;
	let stemEnd = snvIndex + SNV_BASE_BUFFER;

	// 3. Loop to grow stem until we go above desired Tm OR we've come up against both primers
	while (true) {
		// 3a. Slice current stem
		const currentStem = targetStrandSeqSnapPrimerRefPoint.slice(
			stemStart,
			stemEnd + 1,
		);
		// 3b. Get the currentVariantStem
		const currentVariantStem =
			targetStrandSeqSnapPrimerRefPoint.slice(
				stemStart,
				snvSiteSnapPrimerRefPoint.index,
			) +
			snvSiteSnapPrimerRefPoint.variantBase +
			targetStrandSeqSnapPrimerRefPoint.slice(
				snvSiteSnapPrimerRefPoint.index + 1,
				stemEnd + 1,
			);
		// 3c. Calculating the loop length. The loop-side strong mismatch is
		// only inserted when the primer 5' base and base left of the stem are
		// complementary in this candidate stem placement.
		const loopLen = getSnapbackLoopLength(
			targetStrandSeqSnapPrimerRefPoint,
			stemStart,
		);

		// 3d. Build mismatch object for wild and variant type, if needed, for stem Tm calculation
		let wildMismatch = null;
		let variantMismatch = null;
		if (!matchesWild) {
			wildMismatch = {
				position: snvIndex - stemStart, // Appropriate position is relative to the start of the stem
				type: snapbackTailBaseAtSNV,
			};
		} else {
			variantMismatch = {
				position: snvIndex - stemStart, // Appropriate position is relative to the start of the stem
				type: snapbackTailBaseAtSNV,
			};
		}

		// 3e. Compute the wild type allele melting temperature of the snapback
		const wildTm = await calculateSnapbackTmWittwer(
			currentStem,
			loopLen,
			wildMismatch,
		);

		// 3f. Update the closest to desired wild type snapback melting temperature and corresponding stem location if applicable
		if (
			!bestWildTm ||
			Math.abs(wildTm - targetSnapMeltTemp) <
				Math.abs(bestWildTm - targetSnapMeltTemp)
		) {
			bestWildTm = wildTm;
			bestStemLoc.start = stemStart;
			bestStemLoc.end = stemEnd;
			// Compute and save the corresponding best variant type snapback temperature
			correspondingVariantStemTm = await calculateSnapbackTmWittwer(
				currentVariantStem,
				loopLen,
				variantMismatch,
			);
		}

		// 3g. Loop temination if wildTm has become larger than the desired snapback melting temperature
		// for the wild type allele
		if (wildTm >= targetSnapMeltTemp) {
			break;
		}

		// 3h. Grow the stem in the appropriate direction (if it can be grown without overlapping a primer location)
		if (
			stemStart > primerLen &&
			(snvIndex - stemStart < stemEnd - snvIndex ||
				!(stemEnd < seqLen - compPrimerLen - 1))
		) {
			// 3hI. Push the start of the stem one nucleotide to the left only if (the stem is not going to overlap with
			// the primer attachment location) AND [(the beginning of the stem is closer to the SNV that the end of
			// the stem) OR (the end of the stem is up against the reverse primers attachment location (in this frame
			// of reference))]
			stemStart -= 1;
		} else if (stemEnd < seqLen - compPrimerLen - 1) {
			// 3hII. Otherwise we push the end of the stem one nucleotide if (the end of the stem is not up against the
			// reverse primer attachment location)
			// We should push the start of the stem one nucleotide to the left
			stemEnd += 1;
		} else {
			// 3hIII. If we can do neither, we break out of the loop as the stem as grown as large as it can without
			// interfering with primer attachment locations
			break;
		}
	}

	// 4. Final check if final stem doesn’t meet minimum melting temperature requirement
	if (bestWildTm < MINIMUM_TARGET_SNAPBACK_MELTING_TEMP) {
		throw new Error(
			`Could not meet minimum snapback melting temp of ${MINIMUM_TARGET_SNAPBACK_MELTING_TEMP}°C. Final wildTm = ${bestWildTm.toFixed(
				2,
			)}°C. Please consider moving primers farther out so a larger, more stable snapback stem can be created. `,
		);
	}

	// 5. Return the created stem, with its wild and variant allele snapback melting temperatures.
	return {
		bestStemLoc: bestStemLoc,
		meltingTemps: {
			wildTm: parseFloat(bestWildTm.toFixed(TM_DECIMAL_PLACES)),
			variantTm: parseFloat(
				correspondingVariantStemTm.toFixed(TM_DECIMAL_PLACES),
			),
		},
	};
}

/**
 * Constructs the final snapback primer in the reference frame of the primer
 * receiving the snapback tail.
 *
 * The final sequence is composed of (from the 5' end to the 3' end):
 *   1. A strong mismatch at the stem end, which keeps the complement snapback
 *      from extending on its end. The strong mismatch is therefore a complement
 *      to the strong mismatch on the complement strand at the stem's end
 *   2. The stem region of the snapbacks tail
 *   3. An optional strong inner-loop mismatch to prevent one extra loop-closing pair
 * 	 4. The primer
 *
 * Again, the final string is returned 5' → 3'
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Assumptions
 * ──────────────────────────────────────────────────────────────────────────
 * - All DNA strings and objects are in the frame of the primer that recieves the
 *   snapback tail and are passed as 5' to 3'.
 * - The primerLen and compPrimerLen refer to the primer and limiting primer
 *   lengths on this strand and its complement, respectively.
 * 
 * - The SNV lies within the stem
 * - Additional mismatch bases can be appended without exceeding sequence bounds.
 * 		- This is not tested for as it is possible if the minimum primer length
 * 		  exceeds the number of required inner loop and end-of-stem mismatch
 *        bases.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Parameters, Returns, and Errors
 * ──────────────────────────────────────────────────────────────────────────
 * @typedef {Object} PrimerLensRefPoint
 * @property {number} primerLen      Length of the primer on this strand.
 * @property {number} compPrimerLen  Length of the complementary primer.

 * @typedef {Object} StemLoc
 * @property {number} start          Inclusive start index of the stem.
 * @property {number} end            Inclusive end index of the stem.

 * @typedef {Object} SNVSite
 * @property {number} index          0-based position of the SNV on the sequence.
 * @property {string} variantBase    Variant base at that position ("A", "T", "C", or "G").
 * 
 * 
 * @param {string}         seq            Full sequence of the snapback primer strand (5' → 3')
 * @param {SNVSite}        snv            SNV in this strand’s frame
 * @param {PrimerLensRefPoint} primers    Object with primerLen and compPrimerLen
 * @param {StemLoc}        stem           { start: number, end: number } in this strand
 * @param {string}         tailBaseAtSNV       Complement base on the snapback primers tail at the SNV
 *
 * @returns {string}                      Final snapback primer (3' → 5')
 *
 * @throws {Error}                        If inputs are malformed or out of bounds.
 */
function buildSnapbackAndFinalProducts(seq, snv, primers, stem, tailBaseAtSNV) {
	//──────────────────────────────────────────────────────────────────────────//
	//                            Parameter Checking                            //
	//──────────────────────────────────────────────────────────────────────────//

	// 1. Validate sequence
	if (!isValidDNASequence(seq)) {
		throw new Error(`Invalid sequence: must be uppercase A/T/C/G string.`);
	}

	// 2. Validate SNV object
	if (!isValidSNVObject(snv)) {
		throw new Error(
			`Invalid SNV object: must be { index: int, variantBase: A/T/C/G }. Received: ${JSON.stringify(
				snv,
			)}`,
		);
	}
	if (snv.index >= seq.length) {
		throw new Error(
			`SNV index (${snv.index}) exceeds sequence length (${seq.length}).`,
		);
	}

	// 3. Validate primer lengths
	if (
		typeof primers !== 'object' ||
		primers === null ||
		!Number.isInteger(primers.primerLen) ||
		!Number.isInteger(primers.compPrimerLen) ||
		primers.primerLen < MIN_PRIMER_LEN ||
		primers.compPrimerLen < MIN_PRIMER_LEN
	) {
		throw new Error(
			`primers must be an object with integer primerLen and compPrimerLen ≥ ${MIN_PRIMER_LEN}. Received: ${JSON.stringify(
				primers,
			)}`,
		);
	}
	// check keys
	const allowedPrimerKeys = new Set(['primerLen', 'compPrimerLen']);
	for (const key of Object.keys(primers)) {
		if (!allowedPrimerKeys.has(key)) {
			throw new Error(
				`primerLensSnapPrimerRefPoint contains unexpected key "${key}".`,
			);
		}
	}

	// check lengths do not exceed sequence length
	if (primers.primerLen + primers.compPrimerLen >= seq.length) {
		throw new Error(
			`Primer lengths (${primers.primerLen} + ${primers.compPrimerLen}) cannot equal/exceed sequence length (${seq.length}).`,
		);
	}

	// 4. Validate tail base
	if (
		typeof tailBaseAtSNV !== 'string' ||
		tailBaseAtSNV.length !== 1 ||
		!VALID_BASES.has(tailBaseAtSNV)
	) {
		throw new Error(
			`tailBaseAtSNV must be a single base "A", "T", "C", or "G". Received: "${tailBaseAtSNV}".`,
		);
	}

	// 5. Validate stem location
	if (
		typeof stem !== 'object' ||
		stem === null ||
		!Number.isInteger(stem.start) ||
		!Number.isInteger(stem.end) ||
		stem.start < 0 ||
		stem.end < 0 ||
		stem.start > stem.end ||
		stem.end >= seq.length
	) {
		throw new Error(
			`stem must be { start: int, end: int } with 0 ≤ start ≤ end < seq.length. Received: ${JSON.stringify(
				stem,
			)}`,
		);
	}
	// check keys
	const allowedStemKeys = new Set(['start', 'end']);
	for (const key of Object.keys(stem)) {
		if (!allowedStemKeys.has(key)) {
			throw new Error(`stemLoc contains unexpected key "${key}".`);
		}
	}

	//──────────────────────────────────────────────────────────────────────────//
	//								Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//

	// 0. Aliases and destructuring
	const { primerLen, compPrimerLen } = primers;
	const { start: stemStart, end: stemEnd } = stem;
	const snvIndex = snv.index;

	// 1. Initialize the snapback Primer and the products as descriptive objects
	//	  All of these are 5' to 3'
	let snapback = '';
	const descriptiveUnExtendedSnapbackPrimer = {
		fivePrimerLimSnapExtMismatches: '',
		fivePrimeStem: '',
		fivePrimeInnerLoopMismatches: '',
		forwardPrimer: '',
	};
	const descriptiveExtendedSnapback = {
		fivePrimerLimSnapExtMismatches: '',
		fivePrimeStem: '',
		fivePrimeInnerLoopMismatches: '',
		stuffBetween: '', // Includes most if not all of forward primer, and more (also maybe)
		threePrimeInnerLoopMismatches: '',
		threePrimeStem: '',
		threePrimerLimSnapExtMismatches: '',
		threePrimerRestOfAmplicon: '', // Includes most if not all of reverse complement of reverse primer, and more (also maybe)
		// SNV in 3' stem's reference point
		snvOnThreePrimeStem: {
			indexInThreePrimeStem: -1,
			wildBase: '',
			variantBase: '',
		},
		// SNV in 5' stem's reference point
		snvOnFivePrimeStem: {
			indexInFivePrimeStem: -1,
			tailBaseAtSNV: '',
			matchesWild: false,
			matchesVariant: false,
			compWildBase: '',
			compVariantBase: '',
		},
	};

	const descriptiveExendedLimSnapback = {
		threePrimerLimSnapExtMismatches: '',
		threePrimeStem: '',
		threePrimeInnerLoopMismatches: '',
		stuffBetween: '', // Includes most if not all of rev complement of forward primer, and more (also maybe)
		fivePrimeInnerLoopMismatches: '',
		fivePrimeStem: '',
		fivePrimerLimSnapExtMismatches: '',
		fivePrimerRestOfAmplicon: '', // Includes most if not all of reverse primer, and more (also maybe)
		// SNV in 3' stem's reference point
		snvOnThreePrimeStem: {
			indexInThreePrimeStem: -1,
			wildBase: '',
			variantBase: '',
		},
		// SNV in 5' stem's reference point
		snvOnFivePrimeStem: {
			indexInFivePrimeStem: -1,
			tailBaseAtSNV: '',
			matchesWild: false,
			matchesVariant: false,
			compWildBase: '',
			compVariantBase: '',
		},
	};

	// 2. Create the strong mismatch that prevents extension on the 3' end of the complement
	//    snapback primer
	let fivePrimerLimSnapExtMismatches = '';
	for (
		let i = stemEnd + END_OF_STEM_NUMBER_OF_STRONG_BASE_MISMATCHES_REQUIRED;
		i > stemEnd;
		i--
	) {
		// 2a. Get the base on the strand of the snapback primer
		const base = seq[i];

		// 2b. Get the base on the complementary strand
		const compBase = NUCLEOTIDE_COMPLEMENT[base];

		// 2c. We want to mismatch that base — get a strong mismatch for the complement strand so that
		// 	   the complement snapback primer does not extend on the 3' end
		const mismatchAgainstComp = STRONG_NUCLEOTIDE_MISMATCH[compBase];

		// 2d. We insert the COMPLEMENT of the mismatch into the snapback strand
		const mismatchBaseToInsert = NUCLEOTIDE_COMPLEMENT[mismatchAgainstComp];

		// 2e. Insert mismatch (5' → 3')
		fivePrimerLimSnapExtMismatches += mismatchBaseToInsert;
	}

	// 3. Create the 5' end stem region to the snapback primer.
	//    We are adding the reverse complement of the stem region on the sequence strand in the snapback
	//    primer's reference point.
	//    At the SNV site, insert the base chosen earlier (tailBaseAtSNV).
	let fivePrimeStem = '';
	for (let i = stemEnd; i >= stemStart; i--) {
		// If this is the SNV position, insert the selected tail base
		const baseToInsert =
			i === snvIndex ? tailBaseAtSNV : NUCLEOTIDE_COMPLEMENT[seq[i]];

		fivePrimeStem += baseToInsert;
	}

	const innerLoopMismatchCount = getInnerLoopMismatchCount(seq, stemStart);

	// 4. Create the optional strong mismatch in the inner-loop region (before the stem) (5'->3').
	//    This prevents one extra loop-closing pair only when the primer 5' base
	//    and the base immediately left of the stem would otherwise complement.
	let fivePrimeInnerLoopMismatches = '';
	for (
		let i = stemStart - 1;
		i >= stemStart - innerLoopMismatchCount;
		i--
	) {
		const base = seq[i];
		const mismatch = STRONG_NUCLEOTIDE_MISMATCH[base];
		fivePrimeInnerLoopMismatches += mismatch;
	}

	// 5. Grab the forward primer sequence
	const forwardPrimer = seq.slice(0, primerLen);

	// 6. Assemble the unextended snapback primer and fill its descriptive JSON.
	snapback += fivePrimerLimSnapExtMismatches;
	snapback += fivePrimeStem;
	snapback += fivePrimeInnerLoopMismatches;
	snapback += forwardPrimer;

	descriptiveUnExtendedSnapbackPrimer.fivePrimerLimSnapExtMismatches =
		fivePrimerLimSnapExtMismatches;
	descriptiveUnExtendedSnapbackPrimer.fivePrimeStem = fivePrimeStem;
	descriptiveUnExtendedSnapbackPrimer.fivePrimeInnerLoopMismatches =
		fivePrimeInnerLoopMismatches;
	descriptiveUnExtendedSnapbackPrimer.forwardPrimer = forwardPrimer;

	// 7) Build the extended snapback descriptive parts directly from `seq`.
	//    Let the stem location bound everything. Immediately before the stem
	//    is the optional 3′ inner-loop mismatch. Immediately after the stem are the
	// 	  complements to strong mismatches that prevent extension on the
	//    reverse-complement snapback. Everything between the end of the forward
	//    primer and the left inner-loop block is "stuffBetween"; everything after
	//    the right-side end-mismatch block is "threePrimerRestOfAmplicon".

	// Indices around the left-side (5′) inner-loop block
	const innerLoopMismatchesStartLoc =
		stemStart - innerLoopMismatchCount; // inclusive

	// stuffBetween: includes the forward primer and every base up to (but not including)
	// the first 5′ inner-loop mismatch base at innerLoopMismatchesStartLoc.
	const stuffBetween = seq.slice(0, innerLoopMismatchesStartLoc);

	// Inner loop mismatches are right next to the stem
	const threePrimeInnerLoopMismatches = seq.slice(
		innerLoopMismatchesStartLoc,
		stemStart,
	);

	const threePrimeStem = seq.slice(stemStart, stemEnd + 1);

	// threePrimerLimSnapExtMismatches are right to the right of the stem
	const threePrimerLimSnapExtMismatches = seq.slice(
		stemEnd + 1,
		stemEnd + 1 + END_OF_STEM_NUMBER_OF_STRONG_BASE_MISMATCHES_REQUIRED,
	);

	const threePrimerRestOfAmplicon = seq.slice(
		stemEnd + 1 + END_OF_STEM_NUMBER_OF_STRONG_BASE_MISMATCHES_REQUIRED,
	);

	// 7a) SNV annotations for the extended products.
	// Canonical SNV index is relative to descriptiveExtendedSnapback.threePrimeStem.
	const wildBaseAtSNV = seq[snvIndex];
	const variantBaseAtSNV = snv.variantBase;

	const indexInThreePrimeStem = snvIndex - stemStart; // 0-based within threePrimeStem
	const indexInFivePrimeStem = stemEnd - snvIndex; // 0-based within fivePrimeStem (reverse order)

	const compWildBaseAtSNV = NUCLEOTIDE_COMPLEMENT[wildBaseAtSNV];
	const compVariantBaseAtSNV = NUCLEOTIDE_COMPLEMENT[variantBaseAtSNV];

	const tailMatchesWild = tailBaseAtSNV === compWildBaseAtSNV;
	const tailMatchesVariant = tailBaseAtSNV === compVariantBaseAtSNV;

	// Attach SNV metadata to the extended snapback descriptor.
	// NOTE: indexInThreePrimeStem is the canonical reference (relative to descriptiveExtendedSnapback.threePrimeStem).
	descriptiveExtendedSnapback.snvOnThreePrimeStem = {
		indexInThreePrimeStem: indexInThreePrimeStem,
		wildBase: wildBaseAtSNV,
		variantBase: variantBaseAtSNV,
	};

	descriptiveExtendedSnapback.snvOnFivePrimeStem = {
		indexInFivePrimeStem: indexInFivePrimeStem,
		tailBaseAtSNV: tailBaseAtSNV,
		matchesWild: tailMatchesWild,
		matchesVariant: tailMatchesVariant,
		compWildBase: compWildBaseAtSNV,
		compVariantBase: compVariantBaseAtSNV,
	};

	// Fill the rest of descriptive object using the established naming
	descriptiveExtendedSnapback.fivePrimerLimSnapExtMismatches =
		fivePrimerLimSnapExtMismatches;
	descriptiveExtendedSnapback.fivePrimeStem = fivePrimeStem;
	descriptiveExtendedSnapback.fivePrimeInnerLoopMismatches =
		fivePrimeInnerLoopMismatches;
	descriptiveExtendedSnapback.stuffBetween = stuffBetween;
	descriptiveExtendedSnapback.threePrimeInnerLoopMismatches =
		threePrimeInnerLoopMismatches;
	descriptiveExtendedSnapback.threePrimeStem = threePrimeStem;
	descriptiveExtendedSnapback.threePrimerLimSnapExtMismatches =
		threePrimerLimSnapExtMismatches;
	descriptiveExtendedSnapback.threePrimerRestOfAmplicon =
		threePrimerRestOfAmplicon;

	// 8) Extended limiting snapback descriptors by symmetry (reverse-complement)
	const rc = (s) => (s ? reverseComplement(s) : '');

	descriptiveExendedLimSnapback.threePrimerLimSnapExtMismatches = rc(
		fivePrimerLimSnapExtMismatches,
	);
	descriptiveExendedLimSnapback.threePrimeStem = rc(fivePrimeStem);
	descriptiveExendedLimSnapback.threePrimeInnerLoopMismatches = rc(
		fivePrimeInnerLoopMismatches,
	);
	descriptiveExendedLimSnapback.stuffBetween = rc(stuffBetween);
	descriptiveExendedLimSnapback.fivePrimeInnerLoopMismatches = rc(
		threePrimeInnerLoopMismatches,
	);
	descriptiveExendedLimSnapback.fivePrimeStem = rc(threePrimeStem);
	descriptiveExendedLimSnapback.fivePrimerLimSnapExtMismatches = rc(
		threePrimerLimSnapExtMismatches,
	);
	descriptiveExendedLimSnapback.fivePrimerRestOfAmplicon = rc(
		threePrimerRestOfAmplicon,
	);
	// 8a) SNV indices/bases on the extended *limiting* snapback
	descriptiveExendedLimSnapback.snvOnThreePrimeStem = {
		indexInThreePrimeStem: indexInThreePrimeStem,
		wildBase: wildBaseAtSNV,
		variantBase: variantBaseAtSNV,
	};

	descriptiveExendedLimSnapback.snvOnFivePrimeStem = {
		indexInFivePrimeStem: threePrimeStem.length - 1 - indexInThreePrimeStem,
		tailBaseAtSNV: '',
		matchesWild: false,
		matchesVariant: false,
		compWildBase: NUCLEOTIDE_COMPLEMENT[wildBaseAtSNV],
		compVariantBase: NUCLEOTIDE_COMPLEMENT[variantBaseAtSNV],
	};

	// 9) Return only requested artifacts
	return {
		snapback,
		descriptiveUnExtendedSnapbackPrimer,
		descriptiveExtendedSnapback,
		descriptiveExendedLimSnapback,
	};
}

/*****************************************************************************************/
/************************************ Helper Function ************************************/
/*****************************************************************************************/

function shouldInsertInnerLoopMismatch(seq, stemStart) {
	if (!isValidDNASequence(seq)) {
		throw new Error(`Invalid sequence: must be uppercase A/T/C/G string.`);
	}
	if (
		typeof stemStart !== 'number' ||
		!Number.isInteger(stemStart) ||
		stemStart <= 0 ||
		stemStart > seq.length
	) {
		throw new Error(
			`stemStart must be an integer in [1, ${seq.length}]. Received: ${stemStart}.`,
		);
	}

	const primerFivePrimeBase = seq[0];
	const firstFivePrimeOverhangBase = seq[stemStart - 1];

	return (
		NUCLEOTIDE_COMPLEMENT[primerFivePrimeBase] ===
		firstFivePrimeOverhangBase
	);
}

function getInnerLoopMismatchCount(seq, stemStart) {
	return shouldInsertInnerLoopMismatch(seq, stemStart)
		? INNER_LOOP_NUMBER_OF_STRONG_BASE_MISMATCHES_REQUIRED
		: 0;
}

function getSnapbackLoopLength(seq, stemStart) {
	return stemStart + getInnerLoopMismatchCount(seq, stemStart);
}

/**
 * Determines whether a single-nucleotide variant (SNV) lies too close to the
 * primers for a snapback stem to include the required `SNV_BASE_BUFFER`
 * matched bases on either side.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Assumptions
 * ──────────────────────────────────────────────────────────────────────────
 * - The sequence given by `seqLen` represents the target strand (5'→3').
 * - `primerLen` is the length of the primer binding to that strand.
 * - `compPrimerLen` is the length of the primer on the complementary strand.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Parameters, Returns, and Errors
 * ──────────────────────────────────────────────────────────────────────────
 * @param {number} snvIndex       0-based index of the SNV on the target strand.
 * @param {number} primerLen      Length of the primer on the target strand.
 * @param {number} compPrimerLen  Length of the primer on the complementary strand.
 * @param {number} seqLen         Full length of the target sequence.
 *
 * @returns {boolean}             `true`  → SNV is within `SNV_BASE_BUFFER`
 *                                          of a primer (too close)
 *                                `false` → SNV is safely distant.
 *
 * @throws {Error}                If any argument is missing, non-numeric,
 *                                negative, or out of bounds.
 */
function snvTooCloseToPrimer(snvIndex, primerLen, compPrimerLen, seqLen) {
	//──────────────────────────────────────────────────────────────────────────//
	// Parameter Checking                                                      //
	//──────────────────────────────────────────────────────────────────────────//

	// 1. Validating all inputs
	for (const [name, val] of [
		['snvIndex', snvIndex],
		['primerLen', primerLen],
		['compPrimerLen', compPrimerLen],
		['seqLen', seqLen],
	]) {
		if (typeof val !== 'number' || !Number.isFinite(val)) {
			throw new Error(`${name} must be a finite number.`);
		}
		if (!Number.isInteger(val)) {
			throw new Error(`${name} must be an integer.`);
		}
		if (val < 0) {
			throw new Error(`${name} must be non-negative.`);
		}
	}

	if (snvIndex >= seqLen) {
		throw new Error(
			`snvIndex (${snvIndex}) is out of bounds for sequence length ${seqLen}.`,
		);
	}

	if (primerLen + compPrimerLen >= seqLen) {
		throw new Error(
			`primerLen (${primerLen}) + compPrimerLen (${compPrimerLen}) ` +
				`cannot equal or exceed seqLen (${seqLen}).`,
		);
	}

	//──────────────────────────────────────────────────────────────────────────//
	// Function Logic                                                          //
	//──────────────────────────────────────────────────────────────────────────//

	// 1. Calculate the allowable SNV range
	const lowerBound = primerLen + SNV_BASE_BUFFER;
	const upperBound = seqLen - compPrimerLen - SNV_BASE_BUFFER - 1;

	// 2. Return whether the SNV violates either bound
	return snvIndex < lowerBound || snvIndex > upperBound;
}

/**
 * Constructs the mmseq string needed by the Tm service so it sees the intended
 * mismatch in the final double-stranded structure.
 *
 * For example, if mismatch.type = 'G' and the sequence has an A at that location,
 * that means you want an A↔G mismatch in the final pairing. The Tm service expects
 * to see the difference as:
 *   seq=... 'A' ...
 *   mmseq=... 'C' ... (the complement of 'G') at that same position.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Assumptions
 * ──────────────────────────────────────────────────────────────────────────
 * - `seq` is a valid, non-empty, uppercase DNA string (A/T/C/G).
 * - `mismatch` passes `isValidMismatchObject`, meaning:
 *     • `mismatch.position` is a non-negative integer < `seq.length`.
 *     • `mismatch.type` is one of "A", "T", "C", or "G".
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Parameters, Returns, and Errors
 * ──────────────────────────────────────────────────────────────────────────
 * @param {string}   seq       Reference (matched) sequence, 5'→3'.
 * @param {Mismatch} mismatch  { position: number, type: string }
 *
 * @returns {string}           The `mmseq` string to supply to the API.
 *
 * @throws {Error}             If inputs are invalid or out of bounds.
 */
function buildMismatchSequenceForAPI(seq, mismatch) {
	//──────────────────────────────────────────────────────────────────────────//
	//     Parameter Checking                                                   //
	//──────────────────────────────────────────────────────────────────────────//

	// 1. Validate the sequence
	if (!isValidDNASequence(seq)) {
		throw new Error(
			`Invalid DNA sequence: "${seq}". Must contain only A, T, C, or G.`,
		);
	}

	// 2. Validate the mismatch object
	if (!isValidMismatchObject(mismatch)) {
		throw new Error(`Invalid mismatch object: ${JSON.stringify(mismatch)}`);
	}

	// 3. Ensure the mismatch position is within sequence bounds
	if (mismatch.position >= seq.length) {
		throw new Error(
			`Mismatch position (${mismatch.position}) exceeds sequence length ${seq.length}.`,
		);
	}

	//──────────────────────────────────────────────────────────────────────────//
	//     Function Logic                                                       //
	//──────────────────────────────────────────────────────────────────────────//

	// 1. Get the complement of the intended mismatch base
	const compBase = NUCLEOTIDE_COMPLEMENT[mismatch.type];

	// 2. Replace that base in `seq` to produce the mmseq string
	return (
		seq.slice(0, mismatch.position) +
		compBase +
		seq.slice(mismatch.position + 1)
	);
}

/**
 * Parses a raw HTML string to extract a melting temperature (Tm).
 *
 * The input string should include either a <tm> or <mmtm> tag.
 *
 * Example input:
 *   <html><body><seq>...</seq><tm>47.27</tm><mmtm>37.54</mmtm></body></html>
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Assumptions
 * ──────────────────────────────────────────────────────────────────────────
 * - The HTML contains only one <tm> or <mmtm> tag.
 * - The <tm> tag is used for wild-type; <mmtm> is for mismatched variant.
 * - The mismatch flag is optional, but if present, must be a boolean.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Parameters and Returns
 * ──────────────────────────────────────────────────────────────────────────
 * @param {string} rawHtml - Raw HTML string returned from the .cgi file.
 * @param {boolean} [mismatch] - If truthy, extract <mmtm>; otherwise extract <tm>.
 *
 * @returns {number|null} - The extracted Tm as a float, or null if invalid or not found.
 *
 * @throws {Error} - If rawHtml is not a string or mismatch is not a boolean.
 */
function parseTmFromResponse(rawHtml, mismatch) {
	//──────────────────────────────────────────────────────────────────────────//
	//							1. Parameter Checking							//
	//──────────────────────────────────────────────────────────────────────────//

	// 1. Validate rawHTML is a string
	if (typeof rawHtml !== 'string') {
		throw new Error(
			`rawHtml must be a string. Received: ${typeof rawHtml}`,
		);
	}

	// 2. Validate mismatch, if it is passed
	if (mismatch !== undefined && typeof mismatch !== 'boolean') {
		throw new Error(
			`mismatch must be a boolean if provided. Received: ${typeof mismatch}`,
		);
	}

	//──────────────────────────────────────────────────────────────────────────//
	//							2. Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//

	try {
		// 1. Parse the text into a DOM (there might be a better way to do this with functionality built into Node too)
		const parser = new DOMParser();
		const doc = parser.parseFromString(rawHtml, 'text/html');

		console.log('TM PARSING HTML', rawHtml);

		// 2. Getting the <tm> or <mmtm> element
		var tmElement;
		if (!mismatch) {
			// Look for a <tm> element
			tmElement = doc.querySelector('tm');
		} else {
			tmElement = doc.querySelector('mmtm');
		}

		// 3. Returns null if correct tm is not found
		if (!tmElement) {
			return null;
		}

		// 4. Convert the text inside element to a float
		const tmValue = parseFloat(tmElement.textContent.trim());

		// 5. Return parsed Tm, or null if NaN
		return isNaN(tmValue) ? null : tmValue;
	} catch (err) {
		// 6. Fallback: return null on DOM parsing failure
		console.error('parseTmFromResponse error:', err);
		return null;
	}
}

/**
 * Parses ΔH°, ΔS°, and the salt-correction term from a raw HTML response string
 * returned by the dna-utah Tm API.
 *
 * Expected tags in the HTML body:
 *   <dH> -148000.0 </dH>             			(cal/mol)
 *   <dS> -410.0 </dS>                			(cal/K/mol)
 *   <saltCorrection> -11.32 </saltCorrection>  (cal/K/mol)
 *
 * Units:
 *   - dH is converted from cal/mol to kcal/mol before returning.
 *   - dS is returned in cal/(K·mol).
 *   - saltCorrection is returned in cal/K/mol
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Assumptions
 * ──────────────────────────────────────────────────────────────────────────
 * - `rawHtml` form from API keeps tag on top level and the same.
 * - Tags appear at most once; numeric payloads may include whitespace.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Parameters, Returns, and Errors
 * ──────────────────────────────────────────────────────────────────────────
 * @param   {string} rawHtml   Raw HTML response as text.
 *
 * @returns {{ dH:number, dS:number, saltCorrection:number }}
 *          dH in kcal/mol, dS in cal/K/mol, saltCorrection in cal/K/mol
 *
 * @throws  {Error}
 *   - If `rawHtml` is not a string.
 *   - If any of the tags <dH>, <dS>, or <saltCorrection> is missing
 *     or cannot be parsed into a finite number.
 */
function parseThermoParamsFromResponse(rawHtml) {
	//──────────────────────────────────────────────────────────────────────//
	// Parameter Checking                                                   //
	//──────────────────────────────────────────────────────────────────────//
	if (typeof rawHtml !== 'string' || rawHtml.length === 0) {
		throw new Error(
			`rawHtml must be a non-empty string. Received: ${typeof rawHtml}`,
		);
	}

	//──────────────────────────────────────────────────────────────────────//
	// Function Logic                                                       //
	//──────────────────────────────────────────────────────────────────────//
	try {
		// 1. Normalize the input:
		//    - unwrap JSON { contents: "<html>...</html>" }
		//    - drop any preface before first '<'
		//    - fix JSON-escaped closing tags: <\/tag> -> </tag>
		console.log('PARSING THERMO PARAMS RESPONSE HTML', rawHtml);

		let html = rawHtml;
		const trimmed = html.trim();
		if (trimmed.startsWith('{') || trimmed.startsWith('[')) {
			try {
				const obj = JSON.parse(trimmed);
				if (obj && typeof obj.contents === 'string') {
					html = obj.contents;
				}
			} catch {}
		}
		const firstTag = html.indexOf('<');
		if (firstTag > 0) html = html.slice(firstTag);
		html = html.replace(/<\\\//g, '</');

		// 2. Parse HTML response
		const parser = new DOMParser();
		const doc = parser.parseFromString(html, 'text/html');

		// 3. Get each tag and check existence and values
		const dHNode = doc.querySelector('dH');
		const dSNode = doc.querySelector('dS');
		const saltNode = doc.querySelector('saltCorrection');

		if (
			!dHNode ||
			dHNode.textContent == null ||
			!Number.isFinite(Number(dHNode.textContent.trim()))
		) {
			throw new Error('dH not found or unparsable in server response.');
		}
		if (
			!dSNode ||
			dSNode.textContent == null ||
			!Number.isFinite(Number(dSNode.textContent.trim()))
		) {
			throw new Error('dS not found or unparsable in server response.');
		}
		if (
			!saltNode ||
			saltNode.textContent == null ||
			!Number.isFinite(Number(saltNode.textContent.trim()))
		) {
			throw new Error(
				'saltCorrection not found or unparsable in server response.',
			);
		}

		// 4. Parse each tag into a float
		const dH_cal = parseFloat(dHNode.textContent.trim()); // cal/mol
		const dS_cal_per_K = parseFloat(dSNode.textContent.trim()); // cal/K/mol
		const saltCorr_C = parseFloat(saltNode.textContent.trim()); // cal/K/mol

		// 5. Convert delta H to kcal/mol
		const dH_kcal = dH_cal / 1000;

		// 6. Return the thermo parameters
		return { dH: dH_kcal, dS: dS_cal_per_K, saltCorrection: saltCorr_C };
	} catch (err) {
		// Pass along error
		console.error('parseThermoParamsFromResponse error:', err);
		throw err instanceof Error ? err : new Error(String(err));
	}
}

/**
 * Estimates the melting temperature (Tm) of a snapback structure using:
 *
 *     Tm = -5.25 * ln(loopLen) + 0.837 * stemTm + 32.9
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Assumptions
 * ──────────────────────────────────────────────────────────────────────────
 * - `stemSeq` is a valid DNA sequence (uppercase A/C/G/T), 5'→3' direction.
 * - `loopLen` is a number ≥ MIN_LOOP_LEN.
 * - `mismatch` is optional, and if provided:
 *     - Must be an object with shape { position: number, type: string }
 *     - mismatch.position ∈ [SNV_BASE_BUFFER, stemSeq.length - SNV_BASE_BUFFER - 1]
 *     - mismatch.type ∈ { A, T, C, G }
 * - The returned Tm is in degrees Celsius.
 * - getOligoTm(stemSeq, mismatch) returns a numeric Tm.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Parameters and Returns
 * ──────────────────────────────────────────────────────────────────────────
 * @typedef {Object} Mismatch
 * 		@property {number} position   Zero-based index within `seq`
 * 		@property {string} type       Intended base on the opposite strand (A/T/C/G)
 *
 *
 *
 * @param   {string}    stemSeq     DNA stem sequence (5'->3')
 * 									(could be either strand just as long as its
 * 									5'->3')
 * @param   {number}    loopLen     Loop length in nucleotides
 * @param   {Mismatch} [mismatch]   Optional mismatch object:
 *                                  { position: number, type: string }
 *
 * @returns {Promise<number>}       Estimated melting temperature (Tm)
 *
 * @throws  {Error}                 If parameters are invalid
 */
async function calculateSnapbackTmWittwer(stemSeq, loopLen, mismatch) {
	//──────────────────────────────────────────────────────────────────────────//
	//							Parameter Checking								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1. Validate DNA sequence
	if (!isValidDNASequence(stemSeq)) {
		throw new Error(`Invalid DNA sequence: "${stemSeq}"`);
	}

	// 2. Validate loopLen
	if (
		typeof loopLen !== 'number' ||
		!Number.isFinite(loopLen) ||
		loopLen < MIN_LOOP_LEN
	) {
		throw new Error(
			`loopLen must be a finite number ≥ ${MIN_LOOP_LEN}. Received: ${loopLen}`,
		);
	}

	// 3. Validate mismatch object if provided
	if (mismatch !== undefined && mismatch !== null) {
		// 3.1 Validate shape and content
		if (!isValidMismatchObject(mismatch)) {
			throw new Error(
				`Invalid mismatch object: ${JSON.stringify(mismatch)}`,
			);
		}

		// 3.2 Validate mismatch.position bounds
		const min = SNV_BASE_BUFFER;
		const max = stemSeq.length - SNV_BASE_BUFFER - 1;
		if (mismatch.position < min || mismatch.position > max) {
			throw new Error(
				`Mismatch.position (${mismatch.position}) must be between ${min} and ${max} (stem length: ${stemSeq.length})`,
			);
		}
	}

	//──────────────────────────────────────────────────────────────────────────//
	//								Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1. Calculate stem Tm from external method
	const stemTm = await getOligoTm(stemSeq, mismatch ?? undefined);

	// 2. Apply snapback Tm formula
	const tm = -5.25 * Math.log10(loopLen) + 0.837 * stemTm + 32.9;

	// 3. Round result
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
				`Must be a non-empty uppercase string containing only characters A, T, C, and/or G.`,
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
				`Must be a non-empty uppercase string containing only characters A, T, C, and/or G.`,
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
 * Determines whether `seqStrand` is a self-complimentary
 *
 * ─────────────────────────────────────────────────────────────────────────────
 * Assumptions
 * ─────────────────────────────────────────────────────────────────────────────
 * - A valid DNA sequence must be uppercase and contain at least one nucleotide.
 *
 * ─────────────────────────────────────────────────────────────────────────────
 * Parameters, Returns, and Errors
 * ─────────────────────────────────────────────────────────────────────────────
 * @param   {string}  seqStrand   DNA sequence
 *
 * @returns {boolean}             true  → sequence is self-complimentary
 *                                false → sequence is not self-complimentary
 */
function isSelfComplimentary(seqStrand) {
	//──────────────────────────────────────────────────────────────────────────//
	//							Parameter Checking								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1. seqStrand
	if (!isValidDNASequence(seqStrand)) {
		throw new Error(
			`Invalid DNA sequence: ${seqStrand}. ` +
				`Must be a non-empty uppercase string containing only characters A, T, C, and/or G.`,
		);
	}
	//──────────────────────────────────────────────────────────────────────//
	//                           Function Logic                             //
	//──────────────────────────────────────────────────────────────────────//
	// 1. If sequence is odd-length, it can't be self-complimentary
	if (seqStrand.length % 2 === 1) return false;

	// 2. Checks self-complementation base by base
	for (let i = 0, j = seqStrand.length - 1; i < j; i++, j--) {
		if (seqStrand[i] !== NUCLEOTIDE_COMPLEMENT[seqStrand[j]]) {
			return false;
		}
	}

	// 3. If we passed all other checks, the sequence is self-complimentary
	return true;
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
 * Determines whether a given object is a valid mismatch specification.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Assumptions
 * ──────────────────────────────────────────────────────────────────────────
 * - Mismatch must be a plain object (not null or an array).
 * - Must contain exactly two keys: `position` and `type`.
 * - `position` must be a non-negative integer.
 * - `type` must be a valid DNA base: "A", "T", "C", or "G".
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Parameters, Returns, and Errors
 * ──────────────────────────────────────────────────────────────────────────
 * @param   {Object}  mismatch     Object to validate as a mismatch spec.
 *
 * @returns {boolean}              true  → valid mismatch object
 *                                 false → otherwise
 */
function isValidMismatchObject(mismatch) {
	//──────────────────────────────────────────────────────────────────────────//
	//							Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1. Type check: Must be a non-null object and not an array
	if (
		typeof mismatch !== 'object' ||
		mismatch === null ||
		Array.isArray(mismatch)
	) {
		return false;
	}

	// 2. Key check: Must contain exactly 'position' and 'type'
	const expectedKeys = new Set(['position', 'type']);
	const actualKeys = Object.keys(mismatch);
	if (actualKeys.length !== expectedKeys.size) {
		return false;
	}
	for (const key of actualKeys) {
		if (!expectedKeys.has(key)) {
			return false;
		}
	}

	// 3. Validate `position`: must be a non-negative integer
	if (
		typeof mismatch.position !== 'number' ||
		!Number.isInteger(mismatch.position) ||
		mismatch.position < 0
	) {
		return false;
	}

	// 4. Validate `type`: must be a valid base
	if (
		typeof mismatch.type !== 'string' ||
		mismatch.type.length !== 1 ||
		!VALID_BASES.has(mismatch.type)
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
				snvSite,
			)}`,
		);
	}

	// 2. Validate seqLen
	if (
		typeof seqLen !== 'number' ||
		!Number.isInteger(seqLen) ||
		seqLen <= 0
	) {
		throw new Error(
			`seqLen must be a positive integer. Received: ${seqLen}`,
		);
	}

	// 3. Ensure index is within sequence bounds
	if (snvSite.index >= seqLen) {
		throw new Error(
			`snvSite.index (${snvSite.index}) is out of bounds for sequence length ${seqLen}`,
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
				`Must be a non-empty uppercase string containing only characters A, T, C, and/or G..`,
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
	calculateMeltingTempDifferences,
	useForwardPrimer,
	evaluateSnapbackTailMatchingOptions,
	getOligoTm,
	getPrimerTm,
	getThermoParams,
	createStem,
	buildSnapbackAndFinalProducts,
	shouldInsertInnerLoopMismatch,
	getInnerLoopMismatchCount,
	getSnapbackLoopLength,

	// Helper/logic functions
	snvTooCloseToPrimer,
	buildMismatchSequenceForAPI,
	parseTmFromResponse,
	parseThermoParamsFromResponse,
	calculateSnapbackTmWittwer,
	calculateSnapbackTmRochester,
	calculateSnapbackTmSantaLucia,

	// DNA utility functions
	isValidDNASequence,
	isValidSNVObject,
	isValidMismatchObject,
	complementSequence,
	reverseComplement,
	revCompSNV,
	reverseSequence,
	isSelfComplimentary,

	// Constants
	SNV_BASE_BUFFER,
	NUCLEOTIDE_COMPLEMENT,
	MIN_LOOP_LEN,
	MIN_PRIMER_LEN,
	END_OF_STEM_NUMBER_OF_STRONG_BASE_MISMATCHES_REQUIRED,
	MAX_AMPLICON_LEN,
};


// Optional Santa Lucia/Rochester data and helper exports live in their own file.
// Comment this block out if script.js should expose only the primary Wittwer API.
export {
	calculateTm,
	getRochesterHairpinLoopParams,
	getSantaLuciaHicksHairpinParams,
	getDanglingEndParams,
	normalizeDanglingOrientation,
	normalizeNNStep,
	buildTerminalMismatchKey,
	parseTerminalMismatchToken,
	getTerminalMismatchParamsFromToken,
	getTerminalMismatchParams,
	HAIRPIN_LOOP_PARAMETER_ROCHESTER,
	HAIRPIN_LOOP_PARAMETERS_SANTA_LUCIA_HICKS,
	DANGLING_END_PARAMS,
	DANGLING_ORIENTATION,
	TERMINAL_MISMATCH_PARAMS,
} from './optionalTmMethods.js';
