/*
File:           script.js
Description:    Main JavaScript file for interactivity, functionality, and processing.
Author:         Adrian deCola
Relative Path:  uSnapback/src/script.js
*/

import {
	DEFAULT_MAGNESIUM_MM,
	DEFAULT_MONOVALENT_MM,
} from './js/shared/constants.js';
import {
	calculateDuplexThermodynamics,
	calculateSnapbackTmSantaLucia as calculateSnapbackTmSantaLuciaInHouse,
	calculateSnapbackTmWittwer as calculateSnapbackTmWittwerInHouse,
	calculateSnapbackTmWittwerFromStructure,
	calculateTmFromThermodynamics,
	getSantaLuciaHairpinLoopParams,
	normalizeInHouseTmConditions,
} from './js/tm/snapbackTm.js';
import { calculateOwczarzySaltCorrection } from './js/tm/saltCorrection.js';

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
const NO_ADMISSIBLE_STEM_CODE = 'NO_ADMISSIBLE_STEM';
// Chemistry parameters
const T_PARAM = 'SantaLuciaHicks';
// The legacy endpoint does not recognize "owczarzy" and falls back to this mode.
const SALT_CALC_TYPE = 'bpdenominator';
const O_TYPE = 'oligo';
const PRIMER_O_TYPE = 'primer';
// The empirical stem-duplex calculation defaults to both strands at 0.5 µM.
// Callers can override these with tmConditions.concentrationUm and
// tmConditions.limitingConcentrationUm.
const CONC = 0.5;
const LIMITING_CONC = 0.5;
// This gets replaced by build.js
const API_URL = __API_URL__;
const PROXY_URL = __PROXY_URL__;
const USE_PROXY = __USE_PROXY__; // true  |  false  (a real Boolean)
const USE_TOKEN = __USE_TOKEN__; // true | false (real Boolean)
const API_TOKEN = __API_TOKEN__; // string or "" (build-time)

// Retained for backward compatibility only. Rochester is no longer part of the
// production design path; the app now uses the in-house SantaLucia calculator.
const ENABLE_OPTIONAL_TM_METHODS = true;

/*****************************************************************************************/
/************************************ Primary Function ***********************************/
/*****************************************************************************************/

/**
 * Creates a snapback primer by:
 *  1. Independently growing all four primer-orientation × allele-match options.
 *     Every candidate starts with `SNV_BASE_BUFFER` matched bases on each side
 *     of the single nucleotide variant and evaluates every viable stem length.
 *  2. For each option, retaining the stem whose wild-type SantaLucia snapback
 *     Tm is closest to `targetSnapMeltTemp` while remaining at or above 40 °C,
 *     then selecting the option with the largest final wild/variant Tm separation.
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
 * - Snapback Tm values are computed locally with the complete SantaLucia/Hicks
 *   model and stem-only Owczarzy salt correction.
 * - We want to miminize the loop length to the primer on which the snapback tail is on.
 *   If the 5′ base of that primer complements the base immediately left of the stem,
 *   one strong mismatch is inserted at the 3′ end of the snapback tail so the loop
 *   does not zip.
 * - We want to keep the SNV in the middle of the stem or as close to centered as we
 *   can, as we build the snapback stem
 * - Extension can occur on the snapback's complement at the non-loop stem end.
 *   Exactly one strong terminal mismatch pair is used to discourage extension.
 *   on its complement snapback
 * - All sequences are expressed 5'→3' in the frame of the primer that receives the tail.
 * - The optional inner-loop strong mismatch is immediately 5' of the stem (left of the stem).
 * - The strong end-of-stem mismatch used to block extension on the complementary snapback
 *   occur immediately 3' of the stem (right of the stem on this strand).
 * - For the extended descriptive object, stuffBetween includes the forward primer and all
 *   subsequent bases up to (but not including) the first inner-loop mismatch base.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Type definitions
 * ──────────────────────────────────────────────────────────────────────────
 * @typedef {Object} SNVSite
 * @property {number} 				index				0-based position of the SNV on `targetSeqStrand`
 * @property {string}				variantBase			Variant base ('A', 'C', 'G', or 'T')
 *
 * @typedef {Object} TmConditions
 * @property {number} magnesiumMm   Free magnesium concentration in mM
 * @property {number} monovalentMm  Total monovalent-cation concentration in mM
 * @property {number} [concentrationUm]          Empirical reference-strand concentration in µM
 * @property {number} [limitingConcentrationUm]  Empirical partner-strand concentration in µM
 * @property {'log10'} [wittwerLogBase]           Empirical loop logarithm convention; only log10 is supported
 *
 * @typedef {Object} SnapbackMeltingTempDiffs
 * @property {{matchWild:number|null, matchVariant:number|null}} onForwardPrimer  ΔTms (°C) if tail is on the forward primer
 * @property {{matchWild:number|null, matchVariant:number|null}} onReversePrimer  ΔTms (°C) if tail is on the reverse primer
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
 * @property {SnapbackMeltingTempDiffs}		meltingTempDiffs 	Wild/variant ΔTm values for four independently target-optimized designs.
 * @property {Object} optimizedSnapbackOptions Component Tms and stem bounds in both receiving-primer and original-input frames for those four designs.
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
 * @param {number}				targetSnapMeltTemp				Whole-number desired wild-type snapback Tm in the web-app range, 40-80 °C
 * @param {TmConditions} [tmConditions]				Ionic conditions; defaults to 2.2 mM free Mg²⁺ and 13.7 mM total monovalent cations
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
	tmConditions,
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
		!Number.isInteger(targetSnapMeltTemp) ||
		targetSnapMeltTemp < MINIMUM_TARGET_SNAPBACK_MELTING_TEMP ||
		targetSnapMeltTemp > MAXIMUM_TARGET_SNAPBACK_MELTING_TEMP
	) {
		throw new Error(
			`targetSnapMeltTemp must be a whole number from ${MINIMUM_TARGET_SNAPBACK_MELTING_TEMP} to ${MAXIMUM_TARGET_SNAPBACK_MELTING_TEMP} °C.`,
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
	const normalizedTmConditions = normalizeTmConditions(tmConditions);

	//──────────────────────────────────────────────────────────────────────────//
	//								Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//

	// Independently optimize forward/reverse × wild/variant-match choices. This
	// avoids assuming that the largest seven-base seed separation remains best
	// after the complete stem and its terminal contexts have been grown.
	const {
		selected,
		meltingTempDiffs,
		optimizedSnapbackOptions,
	} = await optimizeSnapbackDesignOptions(
		targetSeqStrand,
		snvSite,
		{ primerLen, compPrimerLen },
		targetSnapMeltTemp,
		normalizedTmConditions,
	);
	const tailOnForwardPrimer = selected.context.tailOnForwardPrimer;
	const targetStrandSeqSnapPrimerRefPoint = selected.context.seq;
	const primerLensSnapPrimerRefPoint = selected.context.primerLens;
	const matchesWild = selected.matchesWild;
	const meltingTemps = {
		wildTm: selected.santaLucia.wildTm,
		variantTm: selected.santaLucia.variantTm,
	};
	const {
		snapback,
		descriptiveUnExtendedSnapbackPrimer,
		descriptiveExtendedSnapback,
		descriptiveExendedLimSnapback,
	} = selected.products;
	const snapbackTmSantaLucia = selected.santaLucia;
	const snapbackTmWittwer = calculateSnapbackTmWittwerFromStructure(
		descriptiveExtendedSnapback,
		normalizedTmConditions,
	);
	const snapbackTmRochester = null;

	// Return the selected design plus all independently optimized option scores.
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
		optimizedSnapbackOptions,
		tmConditions: normalizedTmConditions,

		descriptiveUnExtendedSnapbackPrimer,
		descriptiveExtendedSnapback,
		descriptiveExendedLimSnapback,

		// Retained compatibility field; Rochester is outside the two-method
		// production path and is intentionally null.
		snapbackTmRochester,

		// Complete component-level SantaLucia result used for the primary Tms.
		snapbackTmSantaLucia,

		// Wittwer/empirical comparison, calculated fully in-house.
		snapbackTmWittwer,
	};
}

async function calculateOptionalSnapbackTms(
	descriptiveExtendedSnapback,
	tmConditions,
) {
	if (!ENABLE_OPTIONAL_TM_METHODS) {
		return {
			snapbackTmRochester: null,
			snapbackTmSantaLucia: null,
			snapbackTmWittwer: null,
		};
	}

	const snapbackTmSantaLucia = await calculateSnapbackTmSantaLucia(
		descriptiveExtendedSnapback,
		tmConditions,
	);
	const snapbackTmWittwer = calculateSnapbackTmWittwerFromStructure(
		descriptiveExtendedSnapback,
		tmConditions,
	);

	return {
		snapbackTmRochester: null,
		snapbackTmSantaLucia,
		snapbackTmWittwer,
	};
}

async function calculateSnapbackTmRochester(
	descriptiveExtendedSnapback,
	tmConditions,
) {
	const optionalTmMethods = await import('./optionalTmMethods.js');
	return optionalTmMethods.calculateSnapbackTmRochester(
		descriptiveExtendedSnapback,
		getOptionalTmMethodOptions(tmConditions),
	);
}

async function calculateSnapbackTmSantaLucia(
	descriptiveExtendedSnapback,
	tmConditions,
) {
	return calculateSnapbackTmSantaLuciaInHouse(
		descriptiveExtendedSnapback,
		tmConditions,
	);
}

function getOptionalTmMethodOptions(tmConditions) {
	return {
		getThermoParams: (seq, concentration, limitingConc, mismatch) =>
			getThermoParams(
				seq,
				concentration,
				limitingConc,
				mismatch,
				tmConditions,
			),
		conc: CONC,
		limitingConc: LIMITING_CONC,
	};
}

/**
 * Independently grow the four orientation × allele-match designs with the
 * complete SantaLucia model, then choose the final design with the largest
 * unrounded wild/variant separation. Each option first gets its own stem that
 * is closest to the requested wild-type Tm while remaining at or above 40 °C.
 */
async function optimizeSnapbackDesignOptions(
	targetSeqStrand,
	snvSite,
	primerLens,
	targetSnapMeltTemp,
	tmConditions,
) {
	const reverseSeq = reverseComplement(targetSeqStrand);
	const reverseSnv = revCompSNV(snvSite, targetSeqStrand.length);
	const contexts = [
		{
			tailOnForwardPrimer: true,
			seq: targetSeqStrand,
			snv: snvSite,
			primerLens,
		},
		{
			tailOnForwardPrimer: false,
			seq: reverseSeq,
			snv: reverseSnv,
			primerLens: {
				primerLen: primerLens.compPrimerLen,
				compPrimerLen: primerLens.primerLen,
			},
		},
	];

	const candidates = [];
	const failures = [];
	let order = 0;
	for (const context of contexts) {
		for (const matchesWild of [true, false]) {
			const tailBaseAtSNV = NUCLEOTIDE_COMPLEMENT[
				matchesWild
					? context.seq[context.snv.index]
					: context.snv.variantBase
			];
			try {
				const { bestStemLoc } = await createStem(
					context.seq,
					context.snv,
					context.primerLens,
					tailBaseAtSNV,
					matchesWild,
					targetSnapMeltTemp,
					tmConditions,
				);
				const products = buildSnapbackAndFinalProducts(
					context.seq,
					context.snv,
					context.primerLens,
					bestStemLoc,
					tailBaseAtSNV,
				);
				const santaLucia = calculateSnapbackTmSantaLuciaInHouse(
					products.descriptiveExtendedSnapback,
					tmConditions,
				);
				const wildUnrounded = santaLucia.alleles.wild.unroundedTm;
				const variantUnrounded = santaLucia.alleles.variant.unroundedTm;
				candidates.push({
					order,
					context,
					matchesWild,
					tailBaseAtSNV,
					bestStemLoc,
					products,
					santaLucia,
					wildUnrounded,
					variantUnrounded,
					deltaTmUnrounded: Math.abs(wildUnrounded - variantUnrounded),
					targetDistance: Math.abs(
						wildUnrounded - targetSnapMeltTemp,
					),
					stemLength: bestStemLoc.end - bestStemLoc.start + 1,
				});
			} catch (error) {
				if (error?.code !== NO_ADMISSIBLE_STEM_CODE) throw error;
				failures.push({
					order,
					tailOnForwardPrimer: context.tailOnForwardPrimer,
					matchesWild,
					highestWildTm: error.highestWildTm,
				});
			}
			order += 1;
		}
	}

	if (candidates.length === 0) {
		const highestWildTm = Math.max(
			...failures.map((failure) => failure.highestWildTm),
		);
		const error = new Error(
			`No orientation or allele-match option reached ${MINIMUM_TARGET_SNAPBACK_MELTING_TEMP}°C. Highest wildTm = ${highestWildTm.toFixed(2)}°C.`,
		);
		error.code = NO_ADMISSIBLE_STEM_CODE;
		error.highestWildTm = highestWildTm;
		error.candidateFailures = failures;
		throw error;
	}

	const selected = [...candidates].sort((left, right) => {
		if (left.deltaTmUnrounded !== right.deltaTmUnrounded) {
			return right.deltaTmUnrounded - left.deltaTmUnrounded;
		}
		if (left.targetDistance !== right.targetDistance) {
			return left.targetDistance - right.targetDistance;
		}
		if (left.stemLength !== right.stemLength) {
			return left.stemLength - right.stemLength;
		}
		// Preserve the legacy exact-tie preference: variant, then reverse.
		return right.order - left.order;
	})[0];

	const meltingTempDiffs = {
		onForwardPrimer: { matchWild: null, matchVariant: null },
		onReversePrimer: { matchWild: null, matchVariant: null },
	};
	const optimizedSnapbackOptions = {
		onForwardPrimer: { matchWild: null, matchVariant: null },
		onReversePrimer: { matchWild: null, matchVariant: null },
	};
	for (const candidate of candidates) {
		const side = candidate.context.tailOnForwardPrimer
			? 'onForwardPrimer'
			: 'onReversePrimer';
		const match = candidate.matchesWild ? 'matchWild' : 'matchVariant';
		const deltaTm = roundToTmPrecision(candidate.deltaTmUnrounded);
		meltingTempDiffs[side][match] = deltaTm;
		const stemStartOnInputStrand = candidate.context.tailOnForwardPrimer
			? candidate.bestStemLoc.start
			: targetSeqStrand.length - 1 - candidate.bestStemLoc.end;
		const stemEndOnInputStrand = candidate.context.tailOnForwardPrimer
			? candidate.bestStemLoc.end
			: targetSeqStrand.length - 1 - candidate.bestStemLoc.start;
		optimizedSnapbackOptions[side][match] = {
			stemStartInPrimerFrame: candidate.bestStemLoc.start,
			stemEndInPrimerFrame: candidate.bestStemLoc.end,
			stemStartOnInputStrand,
			stemEndOnInputStrand,
			stemLength: candidate.stemLength,
			wildTm: candidate.santaLucia.wildTm,
			variantTm: candidate.santaLucia.variantTm,
			deltaTm,
		};
	}

	return {
		selected,
		meltingTempDiffs,
		optimizedSnapbackOptions,
		failures,
	};
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
 * - The complete in-house SantaLucia/Hicks snapback model, with Owczarzy salt
 *   correction applied to the stem, is used to compute all four Tm differences.
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
async function calculateMeltingTempDifferencesLegacy(
	targetStrandSeqSnapPrimerRefPoint,
	snvSiteSnapPrimerRefPoint,
	bestStemLoc,
	tailOnForwardPrimer,
	tmConditions,
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
		calculateSnapbackTmWittwer(
			stemSeqWildSamePrimer,
			loopLenSamePrimer,
			undefined,
			tmConditions,
		), // 0: same-wild baseline
		calculateSnapbackTmWittwer(
			stemSeqWildSamePrimer,
			loopLenSamePrimer,
			variantTailSamePrimerMismatch,
			tmConditions,
		), // 1: same-wild + variant tail

		// Same primer, variant stem
		calculateSnapbackTmWittwer(
			stemSeqVariantSamePrimer,
			loopLenSamePrimer,
			undefined,
			tmConditions,
		), // 2: same-variant baseline
		calculateSnapbackTmWittwer(
			stemSeqVariantSamePrimer,
			loopLenSamePrimer,
			wildTailSamePrimerMismatch,
			tmConditions,
		), // 3: same-variant + wild tail

		// Reverse primer, wild stem
		calculateSnapbackTmWittwer(
			stemSeqWildRevPrimer,
			loopLenRevPrimer,
			undefined,
			tmConditions,
		), // 4: rev-wild baseline
		calculateSnapbackTmWittwer(
			stemSeqWildRevPrimer,
			loopLenRevPrimer,
			variantTailRevPrimerMismatch,
			tmConditions,
		), // 5: rev-wild + variant tail

		// Reverse primer, variant stem
		calculateSnapbackTmWittwer(
			stemSeqVariantRevPrimer,
			loopLenRevPrimer,
			undefined,
			tmConditions,
		), // 6: rev-variant baseline
		calculateSnapbackTmWittwer(
			stemSeqVariantRevPrimer,
			loopLenRevPrimer,
			wildTailRevPrimerMismatch,
			tmConditions,
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
 * Calculate the four displayed option differences with the complete in-house
 * SantaLucia snapback model. Each option is rebuilt as a real snapback so its
 * own loop length, natural/engineered loop-end mismatch, extension-blocking
 * mismatch, internal mismatch, and Owczarzy correction are all represented.
 */
async function calculateMeltingTempDifferences(
	targetStrandSeqSnapPrimerRefPoint,
	snvSiteSnapPrimerRefPoint,
	bestStemLoc,
	tailOnForwardPrimer,
	tmConditions,
	primerLensSnapPrimerRefPoint = undefined,
) {
	if (!isValidDNASequence(targetStrandSeqSnapPrimerRefPoint)) {
		throw new Error('A valid uppercase target DNA sequence is required.');
	}
	if (!isValidSNVObject(snvSiteSnapPrimerRefPoint)) {
		throw new Error('A valid SNV object is required.');
	}
	if (
		!bestStemLoc ||
		!Number.isInteger(bestStemLoc.start) ||
		!Number.isInteger(bestStemLoc.end) ||
		bestStemLoc.start < 0 ||
		bestStemLoc.end >= targetStrandSeqSnapPrimerRefPoint.length ||
		bestStemLoc.start > bestStemLoc.end ||
		snvSiteSnapPrimerRefPoint.index < bestStemLoc.start ||
		snvSiteSnapPrimerRefPoint.index > bestStemLoc.end
	) {
		throw new Error('bestStemLoc must be valid and contain the SNV.');
	}
	if (typeof tailOnForwardPrimer !== 'boolean') {
		throw new Error('tailOnForwardPrimer must be a boolean.');
	}

	const fallbackPrimerLength = Math.min(
		MIN_PRIMER_LEN,
		bestStemLoc.start,
	);
	const primers = primerLensSnapPrimerRefPoint ?? {
		primerLen: fallbackPrimerLength,
		compPrimerLen: fallbackPrimerLength,
	};
	const reverseSeq = reverseComplement(targetStrandSeqSnapPrimerRefPoint);
	const reverseSnv = revCompSNV(
		snvSiteSnapPrimerRefPoint,
		targetStrandSeqSnapPrimerRefPoint.length,
	);
	const reverseStem = {
		start:
			targetStrandSeqSnapPrimerRefPoint.length - bestStemLoc.end - 1,
		end:
			targetStrandSeqSnapPrimerRefPoint.length - bestStemLoc.start - 1,
	};
	const reversePrimers = {
		primerLen: primers.compPrimerLen,
		compPrimerLen: primers.primerLen,
	};

	const evaluate = (seq, snv, primerLengths, stem, tailBaseAtSNV) => {
		const { descriptiveExtendedSnapback } = buildSnapbackAndFinalProducts(
			seq,
			snv,
			primerLengths,
			stem,
			tailBaseAtSNV,
		);
		return calculateSnapbackTmSantaLuciaInHouse(
			descriptiveExtendedSnapback,
			tmConditions,
		);
	};

	const sameWild = evaluate(
		targetStrandSeqSnapPrimerRefPoint,
		snvSiteSnapPrimerRefPoint,
		primers,
		bestStemLoc,
		NUCLEOTIDE_COMPLEMENT[
			targetStrandSeqSnapPrimerRefPoint[snvSiteSnapPrimerRefPoint.index]
		],
	);
	const sameVariant = evaluate(
		targetStrandSeqSnapPrimerRefPoint,
		snvSiteSnapPrimerRefPoint,
		primers,
		bestStemLoc,
		NUCLEOTIDE_COMPLEMENT[snvSiteSnapPrimerRefPoint.variantBase],
	);
	const reverseWild = evaluate(
		reverseSeq,
		reverseSnv,
		reversePrimers,
		reverseStem,
		NUCLEOTIDE_COMPLEMENT[reverseSeq[reverseSnv.index]],
	);
	const reverseVariant = evaluate(
		reverseSeq,
		reverseSnv,
		reversePrimers,
		reverseStem,
		NUCLEOTIDE_COMPLEMENT[reverseSnv.variantBase],
	);

	const same = {
		matchWild: roundToTmPrecision(
			Math.abs(
				sameWild.alleles.wild.unroundedTm -
					sameWild.alleles.variant.unroundedTm,
			),
		),
		matchVariant: roundToTmPrecision(
			Math.abs(
				sameVariant.alleles.wild.unroundedTm -
					sameVariant.alleles.variant.unroundedTm,
			),
		),
	};
	const opposite = {
		matchWild: roundToTmPrecision(
			Math.abs(
				reverseWild.alleles.wild.unroundedTm -
					reverseWild.alleles.variant.unroundedTm,
			),
		),
		matchVariant: roundToTmPrecision(
			Math.abs(
				reverseVariant.alleles.wild.unroundedTm -
					reverseVariant.alleles.variant.unroundedTm,
			),
		),
	};

	return tailOnForwardPrimer
		? { onForwardPrimer: same, onReversePrimer: opposite }
		: { onForwardPrimer: opposite, onReversePrimer: same };
}

function roundToTmPrecision(value) {
	return Number(value.toFixed(TM_DECIMAL_PLACES));
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
async function useForwardPrimerLegacy(targetSeqStrand, snvSite, tmConditions) {
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
			tmConditions,
		);

	// 8) Evaluate Tm differences for snapback tail on complementary strand
	const tailOnReversePrimerScenario =
		await evaluateSnapbackTailMatchingOptions(
			compInitStem,
			mismatchPos,
			revCompSnvSite.variantBase,
			tmConditions,
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
 * Select the primer orientation and allele-matching tail with the complete
 * SantaLucia snapback calculation. The seven-base seed stem retains the
 * required three matched bases on each side of the SNV.
 */
async function useForwardPrimer(
	targetSeqStrand,
	snvSite,
	tmConditions,
	primerLens = undefined,
) {
	if (!isValidDNASequence(targetSeqStrand)) {
		throw new Error('targetSeqStrand must be a valid uppercase DNA sequence.');
	}
	if (!isValidSNVObject(snvSite) || snvSite.index >= targetSeqStrand.length) {
		throw new Error('snvSite must identify a valid position and variant base.');
	}
	if (
		snvSite.index < SNV_BASE_BUFFER ||
		snvSite.index > targetSeqStrand.length - SNV_BASE_BUFFER - 1
	) {
		throw new Error(
			`The SNV needs ${SNV_BASE_BUFFER} matched bases on each side.`,
		);
	}

	// Preserve the older three-argument helper surface. Without primer lengths
	// there is not enough information to construct the two real loops, so this
	// compatibility path compares the two local seed stems only.
	if (primerLens === undefined) {
		const reverseSeq = reverseComplement(targetSeqStrand);
		const reverseSnv = revCompSNV(snvSite, targetSeqStrand.length);
		const forwardStem = targetSeqStrand.slice(
			snvSite.index - SNV_BASE_BUFFER,
			snvSite.index + SNV_BASE_BUFFER + 1,
		);
		const reverseStem = reverseSeq.slice(
			reverseSnv.index - SNV_BASE_BUFFER,
			reverseSnv.index + SNV_BASE_BUFFER + 1,
		);
		const forward = await evaluateSnapbackTailMatchingOptions(
			forwardStem,
			SNV_BASE_BUFFER,
			snvSite.variantBase,
			tmConditions,
		);
		const reverse = await evaluateSnapbackTailMatchingOptions(
			reverseStem,
			SNV_BASE_BUFFER,
			reverseSnv.variantBase,
			tmConditions,
		);
		return forward.bestTmDifference > reverse.bestTmDifference
			? { tailOnForwardPrimer: true, ...forward }
			: { tailOnForwardPrimer: false, ...reverse };
	}

	const forwardPrimers = primerLens;
	if (
		!Number.isInteger(forwardPrimers.primerLen) ||
		!Number.isInteger(forwardPrimers.compPrimerLen)
	) {
		throw new Error('primerLens must contain integer primerLen and compPrimerLen.');
	}

	const forwardStem = {
		start: snvSite.index - SNV_BASE_BUFFER,
		end: snvSite.index + SNV_BASE_BUFFER,
	};
	const reverseSeq = reverseComplement(targetSeqStrand);
	const reverseSnv = revCompSNV(snvSite, targetSeqStrand.length);
	const reverseStem = {
		start: reverseSnv.index - SNV_BASE_BUFFER,
		end: reverseSnv.index + SNV_BASE_BUFFER,
	};
	const reversePrimers = {
		primerLen: forwardPrimers.compPrimerLen,
		compPrimerLen: forwardPrimers.primerLen,
	};

	const evaluateOrientation = (seq, orientedSnv, primers, stem) => {
		const wildTailBase = NUCLEOTIDE_COMPLEMENT[seq[orientedSnv.index]];
		const variantTailBase =
			NUCLEOTIDE_COMPLEMENT[orientedSnv.variantBase];
		const wildStructure = buildSnapbackAndFinalProducts(
			seq,
			orientedSnv,
			primers,
			stem,
			wildTailBase,
		).descriptiveExtendedSnapback;
		const variantStructure = buildSnapbackAndFinalProducts(
			seq,
			orientedSnv,
			primers,
			stem,
			variantTailBase,
		).descriptiveExtendedSnapback;
		const wildResult = calculateSnapbackTmSantaLuciaInHouse(
			wildStructure,
			tmConditions,
		);
		const variantResult = calculateSnapbackTmSantaLuciaInHouse(
			variantStructure,
			tmConditions,
		);
		const wildDifference = Math.abs(
			wildResult.alleles.wild.unroundedTm -
				wildResult.alleles.variant.unroundedTm,
		);
		const variantDifference = Math.abs(
			variantResult.alleles.wild.unroundedTm -
				variantResult.alleles.variant.unroundedTm,
		);
		if (wildDifference > variantDifference) {
			return {
				bestSnapbackTailBaseAtSNV: wildTailBase,
				bestTmDifference: wildDifference,
				snapbackTailMatchesWild: true,
			};
		}
		return {
			bestSnapbackTailBaseAtSNV: variantTailBase,
			bestTmDifference: variantDifference,
			snapbackTailMatchesWild: false,
		};
	};

	const forward = evaluateOrientation(
		targetSeqStrand,
		snvSite,
		forwardPrimers,
		forwardStem,
	);
	const reverse = evaluateOrientation(
		reverseSeq,
		reverseSnv,
		reversePrimers,
		reverseStem,
	);

	if (forward.bestTmDifference > reverse.bestTmDifference) {
		return { tailOnForwardPrimer: true, ...forward };
	}
	return { tailOnForwardPrimer: false, ...reverse };
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
async function evaluateSnapbackTailMatchingOptionsLegacy(
	initStem,
	mismatchPos,
	variantBase,
	tmConditions,
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
	const wildMatchTmPromise = await getOligoTm(
		initStem,
		undefined,
		tmConditions,
	);

	// 2) Get Tm for a stem with the variant allele where the snapback tail matches it
	const variantInitStem =
		initStem.slice(0, mismatchPos) +
		variantBase +
		initStem.slice(mismatchPos + 1);

	const variantMatchTmPromise = await getOligoTm(
		variantInitStem,
		undefined,
		tmConditions,
	);

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
		tmConditions,
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
		tmConditions,
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
 * Backward-compatible seed-stem helper implemented locally with the same
 * SantaLucia/Hicks nearest-neighbour and Owczarzy models as the full calculator.
 * The main web-app path uses complete snapback structures in useForwardPrimer.
 */
async function evaluateSnapbackTailMatchingOptions(
	initStem,
	mismatchPos,
	variantBase,
	tmConditions,
) {
	if (!isValidDNASequence(initStem)) {
		throw new Error('initStem must be a non-empty uppercase DNA sequence.');
	}
	if (
		!Number.isInteger(mismatchPos) ||
		mismatchPos < 1 ||
		mismatchPos >= initStem.length - 1
	) {
		throw new Error('mismatchPos must be an internal stem index.');
	}
	if (!VALID_BASES.has(variantBase) || variantBase === initStem[mismatchPos]) {
		throw new Error('variantBase must be a valid base different from wild type.');
	}

	const normalized = normalizeTmConditions(tmConditions);
	const stemTm = (sequence, mismatch) => {
		const thermo = calculateDuplexThermodynamics(sequence, mismatch);
		const salt = calculateOwczarzySaltCorrection({
			stemSequence: sequence,
			stemDeltaH: thermo.dH,
			magnesiumMm: normalized.magnesiumMm,
			monovalentMm: normalized.monovalentMm,
		});
		return calculateTmFromThermodynamics({
			dH: thermo.dH,
			dS: thermo.dS,
			saltCorrection: salt.saltCorrection,
			concentrationUm: normalized.concentrationUm,
			limitingConcentrationUm: normalized.limitingConcentrationUm,
		});
	};

	const wildBase = initStem[mismatchPos];
	const variantStem =
		initStem.slice(0, mismatchPos) +
		variantBase +
		initStem.slice(mismatchPos + 1);
	const wildMatchTm = stemTm(initStem);
	const variantMatchTm = stemTm(variantStem);
	const wildTailMismatchTm = stemTm(variantStem, {
		position: mismatchPos,
		type: NUCLEOTIDE_COMPLEMENT[wildBase],
	});
	const variantTailMismatchTm = stemTm(initStem, {
		position: mismatchPos,
		type: NUCLEOTIDE_COMPLEMENT[variantBase],
	});
	const wildDifference = Math.abs(wildMatchTm - wildTailMismatchTm);
	const variantDifference = Math.abs(
		variantMatchTm - variantTailMismatchTm,
	);

	if (wildDifference > variantDifference) {
		return {
			bestSnapbackTailBaseAtSNV: NUCLEOTIDE_COMPLEMENT[wildBase],
			bestTmDifference: wildDifference,
			snapbackTailMatchesWild: true,
		};
	}
	return {
		bestSnapbackTailBaseAtSNV: NUCLEOTIDE_COMPLEMENT[variantBase],
		bestTmDifference: variantDifference,
		snapbackTailMatchesWild: false,
	};
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
 * – Ionic conditions default to the values requested by the CTW lab.
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
 * @param  {TmConditions} [tmConditions] Ionic conditions for this calculation
 *
 * @returns {Promise<number>}           Melting temperature (°C) rounded to
 *                                      `TM_DECIMAL_PLACES`
 *
 * @throws {Error}                      If inputs are invalid, the network
 *                                      request fails, or the CGI response
 *                                      cannot be parsed
 */
function normalizeTmConditions(tmConditions) {
	if (
		tmConditions !== undefined &&
		tmConditions !== null &&
		(typeof tmConditions !== 'object' || Array.isArray(tmConditions))
	) {
		throw new Error('tmConditions must be an object when provided.');
	}

	const magnesiumMm =
		tmConditions?.magnesiumMm ?? DEFAULT_MAGNESIUM_MM;
	const monovalentMm =
		tmConditions?.monovalentMm ?? DEFAULT_MONOVALENT_MM;
	const concentrationUm = tmConditions?.concentrationUm ?? CONC;
	const limitingConcentrationUm =
		tmConditions?.limitingConcentrationUm ?? LIMITING_CONC;
	const wittwerLogBase = tmConditions?.wittwerLogBase ?? 'log10';

	for (const [name, value] of [
		['magnesiumMm', magnesiumMm],
		['monovalentMm', monovalentMm],
	]) {
		if (
			typeof value !== 'number' ||
			!Number.isFinite(value) ||
			value < 0
		) {
			throw new Error(`${name} must be a finite, non-negative number.`);
		}
	}
	if (magnesiumMm === 0 && monovalentMm === 0) {
		throw new Error(
			'At least one ionic concentration must be greater than zero for the Owczarzy correction.',
		);
	}
	for (const [name, value] of [
		['concentrationUm', concentrationUm],
		['limitingConcentrationUm', limitingConcentrationUm],
	]) {
		if (
			typeof value !== 'number' ||
			!Number.isFinite(value) ||
			value <= 0
		) {
			throw new Error(`${name} must be a finite, positive number.`);
		}
	}
	if (wittwerLogBase !== 'log10') {
		throw new Error('wittwerLogBase must be "log10".');
	}

	return {
		magnesiumMm,
		monovalentMm,
		concentrationUm,
		limitingConcentrationUm,
		wittwerLogBase,
	};
}

function buildTmRequestParams(
	seq,
	mismatch,
	tmConditions,
	{
		otype = O_TYPE,
		concentration = null,
		limitingConc = null,
	} = {},
) {
	const { magnesiumMm, monovalentMm } =
		normalizeTmConditions(tmConditions);
	// `mg` is the user-supplied free Mg²⁺ value; dNTP is intentionally omitted.
	const params = new URLSearchParams({
		mg: String(magnesiumMm),
		mono: String(monovalentMm),
		seq: seq.toLowerCase(),
		tparam: T_PARAM,
		saltcalctype: SALT_CALC_TYPE,
		otype,
		decimalplaces: String(TM_DECIMAL_PLACES),
	});

	if (concentration != null) {
		params.set('concentration', String(concentration));
	}
	if (limitingConc != null) {
		params.set('limitingconc', String(limitingConc));
	}
	if (USE_TOKEN && API_TOKEN) {
		params.set('token', API_TOKEN);
	}
	if (mismatch) {
		params.set('mmseq', buildMismatchSequenceForAPI(seq, mismatch));
	}

	return params;
}

function buildTmRequestUrl(params) {
	const apiURL = `${API_URL}?${params.toString()}`;
	return USE_PROXY
		? `${PROXY_URL}?url=${encodeURIComponent(apiURL)}`
		: apiURL;
}

async function getOligoTm(seq, mismatch, tmConditions) {
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

	const params = buildTmRequestParams(seq, mismatch, tmConditions, {
		otype: O_TYPE,
		concentration: CONC,
		limitingConc: LIMITING_CONC,
	});
	const finalURL = buildTmRequestUrl(params);

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
 * @param  {TmConditions} [tmConditions] Ionic conditions for this calculation
 * @returns {Promise<number>}           Melting temperature (°C) rounded to
 *                                      `TM_DECIMAL_PLACES`
 * @throws {Error}                      If inputs are invalid, the network
 *                                      request fails, or the CGI response
 *                                      cannot be parsed
 */
async function getPrimerTm(seq, tmConditions) {
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

	const params = buildTmRequestParams(seq, null, tmConditions, {
		otype: PRIMER_O_TYPE,
	});
	const finalURL = buildTmRequestUrl(params);

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
async function getThermoParams(
	seq,
	concentration,
	limitingConc,
	mismatch,
	tmConditions,
) {
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

	const params = buildTmRequestParams(seq, mismatch, tmConditions, {
		concentration: concentration ?? CONC,
		limitingConc: limitingConc ?? LIMITING_CONC,
	});
	const finalURL = buildTmRequestUrl(params);

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
	tmConditions,
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

	// Enlarge the stem right then left while respecting the primer sites, and
	// retain the complete-structure SantaLucia result closest to the requested Tm.

	// 1. Initialize variable to hold the snapback melting temperature for the wild type allele that is closest to the desired
	// snapback melting temperature for the wild type allele. Also initialize a variable for the corresponding melting
	// temperature of the variant snapback melting temperature. Finally, initialize the variable for the corresponding stem locations
	let bestWildTm = null;
	let correspondingVariantStemTm = null;
	let bestStemLoc = { start: null, end: null };
	let highestWildTm = Number.NEGATIVE_INFINITY;

	// 2. Initialize stem region
	const snvIndex = snvSiteSnapPrimerRefPoint.index;
	let stemStart = snvIndex - SNV_BASE_BUFFER;
	let stemEnd = snvIndex + SNV_BASE_BUFFER;

	// 3. Evaluate every viable stem length until both primer boundaries are reached.
	while (true) {
		// 3a. Build this exact candidate and calculate both alleles using the
		// complete intramolecular SantaLucia model. This ensures stem growth sees
		// the loop and both real terminal mismatches, not just an isolated duplex.
		const candidateExtendedSnapback = buildSnapbackAndFinalProducts(
			targetStrandSeqSnapPrimerRefPoint,
			snvSiteSnapPrimerRefPoint,
			primerLensSnapPrimerRefPoint,
			{ start: stemStart, end: stemEnd },
			snapbackTailBaseAtSNV,
		).descriptiveExtendedSnapback;
		const candidateTms = calculateSnapbackTmSantaLuciaInHouse(
			candidateExtendedSnapback,
			tmConditions,
		);
		const wildTm = candidateTms.alleles.wild.unroundedTm;
		highestWildTm = Math.max(highestWildTm, wildTm);

		// 3b. Update the closest admissible wild-type Tm and its corresponding
		// variant Tm. A sub-40 candidate must not hide a viable warmer stem merely
		// because it is numerically closer to a low requested target.
		if (
			wildTm >= MINIMUM_TARGET_SNAPBACK_MELTING_TEMP &&
			(bestWildTm === null ||
				Math.abs(wildTm - targetSnapMeltTemp) <
					Math.abs(bestWildTm - targetSnapMeltTemp))
		) {
			bestWildTm = wildTm;
			bestStemLoc.start = stemStart;
			bestStemLoc.end = stemEnd;
			correspondingVariantStemTm =
				candidateTms.alleles.variant.unroundedTm;
		}

		// 3c. Grow the stem in the appropriate direction (if it can be grown without overlapping a primer location).
		// Evaluate every viable length instead of stopping at the first crossing:
		// sequence-specific terminal-mismatch changes can make the full SantaLucia
		// snapback Tm slightly non-monotonic.
		if (
			stemStart > primerLen &&
			(snvIndex - stemStart < stemEnd - snvIndex ||
				!(stemEnd < seqLen - compPrimerLen - 1))
		) {
			// 3cI. Push the start of the stem one nucleotide to the left only if (the stem is not going to overlap with
			// the primer attachment location) AND [(the beginning of the stem is closer to the SNV that the end of
			// the stem) OR (the end of the stem is up against the reverse primers attachment location (in this frame
			// of reference))]
			stemStart -= 1;
		} else if (stemEnd < seqLen - compPrimerLen - 1) {
			// 3cII. Otherwise we push the end of the stem one nucleotide if (the end of the stem is not up against the
			// reverse primer attachment location)
			// We should push the start of the stem one nucleotide to the left
			stemEnd += 1;
		} else {
			// 3cIII. If we can do neither, we break out of the loop as the stem as grown as large as it can without
			// interfering with primer attachment locations
			break;
		}
	}

	// 4. Final check if no viable stem meets the minimum melting temperature.
	if (bestWildTm === null) {
		const error = new Error(
			`Could not meet minimum snapback melting temp of ${MINIMUM_TARGET_SNAPBACK_MELTING_TEMP}°C. Highest wildTm = ${highestWildTm.toFixed(
				2,
			)}°C. Please consider moving primers farther out so a larger, more stable snapback stem can be created. `,
		);
		error.code = NO_ADMISSIBLE_STEM_CODE;
		error.highestWildTm = highestWildTm;
		throw error;
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
 *     Tm = -5.25 * log10(loopLen) + 0.837 * stemTm + 32.9
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
 * @param   {TmConditions} [tmConditions] Ionic conditions for this calculation
 *
 * @returns {Promise<number>}       Estimated melting temperature (Tm)
 *
 * @throws  {Error}                 If parameters are invalid
 */
async function calculateSnapbackTmWittwer(
	stemSeq,
	loopLen,
	mismatch,
	tmConditions,
) {
	// Complete-structure overload: both public snapback methods can be called as
	// (descriptiveExtendedSnapback, tmConditions). The legacy four-argument
	// stem/loop form remains supported below.
	if (stemSeq && typeof stemSeq === 'object' && !Array.isArray(stemSeq)) {
		return calculateSnapbackTmWittwerFromStructure(
			stemSeq,
			loopLen ?? {},
		);
	}

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

	// Calculate locally with the log10 empirical convention used by the prior
	// uSnapback web app and the reference workbook calculations.
	const normalizedConditions = normalizeTmConditions(tmConditions);
	return calculateSnapbackTmWittwerInHouse(
		stemSeq,
		loopLen,
		mismatch ?? undefined,
		normalizedConditions,
		{ logBase: normalizedConditions.wittwerLogBase },
	);
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
	buildTmRequestParams,
	normalizeTmConditions,
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
	calculateSnapbackTmWittwerFromStructure,
	calculateDuplexThermodynamics,
	calculateTmFromThermodynamics,
	calculateOwczarzySaltCorrection,
	getSantaLuciaHairpinLoopParams,
	normalizeInHouseTmConditions,

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


// Legacy Rochester data and parameter lookup helpers remain available for
// compatibility, but Rochester is not used by the production design path.
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
