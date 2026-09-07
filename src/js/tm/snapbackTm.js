import {
	DEFAULT_MAGNESIUM_MM,
	DEFAULT_MONOVALENT_MM,
} from '../shared/constants.js';
import {
	DNA_COMPLEMENT,
	DUPLEX_INITIATION,
	GAS_CONSTANT_CAL,
	HAIRPIN_LOOP_ANCHORS,
	INTERNAL_MISMATCH_PARAMS,
	MATCHED_NN_PARAMS,
	SYMMETRY_CORRECTION,
	TEMPERATURE_37_K,
	TERMINAL_AT_PENALTY,
	VALID_DNA_BASES,
} from './parameters.js';
import { calculateOwczarzySaltCorrection } from './saltCorrection.js';
import { getTerminalMismatchParams } from '../../optionalTmMethods.js';

export const DEFAULT_STEM_CONCENTRATION_UM = 0.5;
export const DEFAULT_LIMITING_STEM_CONCENTRATION_UM = 0.5;
export const WITTWER_COEFFICIENTS = Object.freeze({
	intercept: 32.9,
	loopLog: -5.25,
	stemTm: 0.837,
	logBase: 'log10',
});

const TM_DECIMAL_PLACES = 2;
const MIN_INTERNAL_MATCHED_FLANK = 3;

function roundTm(value) {
	return Number(value.toFixed(TM_DECIMAL_PLACES));
}

function wittwerLoopLog(loopLength, logBase) {
	if (logBase !== 'log10') {
		throw new Error('The Wittwer/empirical method uses log10 only.');
	}
	return Math.log10(loopLength);
}

/**
 * Combine an already-calculated stem Tm with the Wittwer/empirical regression.
 *
 * This deliberately does not round so callers can retain full precision while
 * selecting a design. The later Hugh-refit workbook also exposes this exact
 * composition step independently of its stem calculator.
 */
export function calculateWittwerTmFromStemTm(
	stemTm,
	loopLength,
	logBase = WITTWER_COEFFICIENTS.logBase,
) {
	if (typeof stemTm !== 'number' || !Number.isFinite(stemTm)) {
		throw new Error('stemTm must be a finite number in degrees C.');
	}
	if (!Number.isInteger(loopLength) || loopLength < 3) {
		throw new Error('loopLength must be an integer of at least 3 bases.');
	}

	return (
		WITTWER_COEFFICIENTS.intercept +
		WITTWER_COEFFICIENTS.loopLog *
			wittwerLoopLog(loopLength, logBase) +
		WITTWER_COEFFICIENTS.stemTm * stemTm
	);
}

function validateDna(name, sequence, minimumLength = 1) {
	if (
		typeof sequence !== 'string' ||
		sequence.length < minimumLength ||
		![...sequence].every((base) => VALID_DNA_BASES.has(base))
	) {
		throw new Error(
			`${name} must be an uppercase A/C/G/T string of length at least ${minimumLength}.`,
		);
	}
}

function reverse(sequence) {
	return [...sequence].reverse().join('');
}

function reverseComplement(sequence) {
	return reverse(sequence)
		.split('')
		.map((base) => DNA_COMPLEMENT[base])
		.join('');
}

function isCanonicalPair(topBase, alignedBottomBase) {
	return DNA_COMPLEMENT[topBase] === alignedBottomBase;
}

function normalizeConditions(conditions = {}) {
	if (!conditions || typeof conditions !== 'object' || Array.isArray(conditions)) {
		throw new Error('Tm conditions must be an object.');
	}
	const magnesiumMm = conditions.magnesiumMm ?? DEFAULT_MAGNESIUM_MM;
	const monovalentMm = conditions.monovalentMm ?? DEFAULT_MONOVALENT_MM;
	const concentrationUm =
		conditions.concentrationUm ?? DEFAULT_STEM_CONCENTRATION_UM;
	const limitingConcentrationUm =
		conditions.limitingConcentrationUm ??
		DEFAULT_LIMITING_STEM_CONCENTRATION_UM;
	const wittwerLogBase =
		conditions.wittwerLogBase ?? WITTWER_COEFFICIENTS.logBase;

	for (const [name, value] of [
		['magnesiumMm', magnesiumMm],
		['monovalentMm', monovalentMm],
	]) {
		if (typeof value !== 'number' || !Number.isFinite(value) || value < 0) {
			throw new Error(`${name} must be a finite, non-negative number.`);
		}
	}
	if (magnesiumMm === 0 && monovalentMm === 0) {
		throw new Error(
			'Owczarzy salt correction requires magnesiumMm or monovalentMm to be greater than zero.',
		);
	}
	for (const [name, value] of [
		['concentrationUm', concentrationUm],
		['limitingConcentrationUm', limitingConcentrationUm],
	]) {
		if (typeof value !== 'number' || !Number.isFinite(value) || value <= 0) {
			throw new Error(`${name} must be a finite, positive number in µM.`);
		}
	}
	if (wittwerLogBase !== 'log10') {
		throw new Error('wittwerLogBase must be "log10".');
	}

	return Object.freeze({
		magnesiumMm,
		monovalentMm,
		concentrationUm,
		limitingConcentrationUm,
		wittwerLogBase,
	});
}

/**
 * Sum SantaLucia/Hicks matched and Allawi/SantaLucia single-mismatch nearest
 * neighbours for two explicit antiparallel stem strands.
 *
 * `topFiveToThree` and `bottomFiveToThree` are both supplied in their own
 * conventional 5'->3' direction. The bottom strand is reversed internally so
 * each character is aligned with its opposing base on the top strand.
 */
export function calculateDuplexThermodynamicsFromStrands(
	topFiveToThree,
	bottomFiveToThree,
	{ includeInitiation = true, includeSymmetry = false } = {},
) {
	validateDna('topFiveToThree', topFiveToThree, 2);
	validateDna('bottomFiveToThree', bottomFiveToThree, 2);
	if (topFiveToThree.length !== bottomFiveToThree.length) {
		throw new Error('The two stem strands must have the same length.');
	}

	const bottomAligned = reverse(bottomFiveToThree);
	const mismatchPositions = [];
	for (let i = 0; i < topFiveToThree.length; i += 1) {
		if (!isCanonicalPair(topFiveToThree[i], bottomAligned[i])) {
			mismatchPositions.push(i);
		}
	}
	if (mismatchPositions.length > 1) {
		throw new Error(
			`Only one internal stem mismatch is supported; found ${mismatchPositions.length}.`,
		);
	}
	if (
		mismatchPositions.length === 1 &&
		(mismatchPositions[0] === 0 ||
			mismatchPositions[0] === topFiveToThree.length - 1)
	) {
		throw new Error(
			'An internal mismatch must have a neighbouring stack on both sides.',
		);
	}

	let dH = includeInitiation ? DUPLEX_INITIATION.dH : 0;
	let dS = includeInitiation ? DUPLEX_INITIATION.dS : 0;
	const stacks = [];

	for (let i = 0; i < topFiveToThree.length - 1; i += 1) {
		const top2 = topFiveToThree.slice(i, i + 2);
		const bottom2 = bottomAligned.slice(i, i + 2);
		const canonical =
			isCanonicalPair(top2[0], bottom2[0]) &&
			isCanonicalPair(top2[1], bottom2[1]);
		const key = `${top2}/${bottom2}`;
		const params = canonical
			? MATCHED_NN_PARAMS[top2]
			: INTERNAL_MISMATCH_PARAMS[key];
		if (!params) {
			throw new Error(`No internal nearest-neighbour parameters for ${key}.`);
		}
		dH += params.dH;
		dS += params.dS;
		stacks.push(
			Object.freeze({
				index: i,
				key,
				kind: canonical ? 'matched' : 'internal-mismatch',
				dH: params.dH,
				dS: params.dS,
			}),
		);
	}

	const terminalAT = [];
	for (const index of [0, topFiveToThree.length - 1]) {
		if (
			isCanonicalPair(topFiveToThree[index], bottomAligned[index]) &&
			(topFiveToThree[index] === 'A' || topFiveToThree[index] === 'T')
		) {
			dH += TERMINAL_AT_PENALTY.dH;
			dS += TERMINAL_AT_PENALTY.dS;
			terminalAT.push(index === 0 ? 'left' : 'right');
		}
	}

	const selfComplementary =
		mismatchPositions.length === 0 &&
		topFiveToThree === reverseComplement(topFiveToThree);
	if (includeSymmetry && selfComplementary) {
		dH += SYMMETRY_CORRECTION.dH;
		dS += SYMMETRY_CORRECTION.dS;
	}

	return Object.freeze({
		dH,
		dS,
		topFiveToThree,
		bottomFiveToThree,
		bottomAligned,
		mismatchPositions: Object.freeze(mismatchPositions),
		stacks: Object.freeze(stacks),
		terminalAT: Object.freeze(terminalAT),
		initiation: includeInitiation ? DUPLEX_INITIATION : null,
		symmetryApplied: includeSymmetry && selfComplementary,
	});
}

/** Build an explicit duplex from the legacy {position, type} mismatch input. */
export function calculateDuplexThermodynamics(
	stemSequence,
	mismatch = undefined,
	options = {},
) {
	validateDna('stemSequence', stemSequence, 2);
	const bottomAligned = [...stemSequence].map((base) => DNA_COMPLEMENT[base]);
	if (mismatch !== undefined && mismatch !== null) {
		if (
			typeof mismatch !== 'object' ||
			!Number.isInteger(mismatch.position) ||
			mismatch.position < 1 ||
			mismatch.position >= stemSequence.length - 1 ||
			!VALID_DNA_BASES.has(mismatch.type)
		) {
			throw new Error(
				'mismatch must be { position: an internal index, type: A/C/G/T }.',
			);
		}
		if (bottomAligned[mismatch.position] === mismatch.type) {
			throw new Error('mismatch.type is complementary and therefore is not a mismatch.');
		}
		bottomAligned[mismatch.position] = mismatch.type;
	}
	return calculateDuplexThermodynamicsFromStrands(
		stemSequence,
		reverse(bottomAligned.join('')),
		options,
	);
}

export function calculateTmFromThermodynamics({
	dH,
	dS,
	saltCorrection = 0,
	concentrationUm,
	limitingConcentrationUm,
}) {
	for (const [name, value] of [
		['dH', dH],
		['dS', dS],
		['saltCorrection', saltCorrection],
	]) {
		if (typeof value !== 'number' || !Number.isFinite(value)) {
			throw new Error(`${name} must be finite.`);
		}
	}

	let concentrationEntropy = 0;
	if (
		concentrationUm !== undefined ||
		limitingConcentrationUm !== undefined
	) {
		if (
			typeof concentrationUm !== 'number' ||
			!Number.isFinite(concentrationUm) ||
			concentrationUm <= 0 ||
			typeof limitingConcentrationUm !== 'number' ||
			!Number.isFinite(limitingConcentrationUm) ||
			limitingConcentrationUm <= 0
		) {
			throw new Error('Both strand concentrations must be positive values in µM.');
		}
		const duplexAtTm =
			Math.min(concentrationUm, limitingConcentrationUm) / 2;
		const freeA = concentrationUm - duplexAtTm;
		const freeB = limitingConcentrationUm - duplexAtTm;
		const concentrationTermM =
			((freeA * freeB) / duplexAtTm) * 1e-6;
		concentrationEntropy = GAS_CONSTANT_CAL * Math.log(concentrationTermM);
	}

	const denominator = dS + saltCorrection + concentrationEntropy;
	const kelvin = (dH * 1000) / denominator;
	if (!Number.isFinite(kelvin) || kelvin <= 0) {
		throw new Error('Thermodynamic inputs produced a non-physical Tm.');
	}
	return kelvin - 273.15;
}

function sumTerminalMismatchParams(terminalMismatches = []) {
	if (!Array.isArray(terminalMismatches)) {
		throw new Error('terminalMismatches must be an array.');
	}
	let dH = 0;
	let dS = 0;
	const entries = terminalMismatches.map(({ top2, bottom2, label }) => {
		const params = getTerminalMismatchParams(top2, bottom2);
		dH += params.dH;
		dS += params.dS;
		return Object.freeze({ label, top2, bottom2, ...params });
	});
	return Object.freeze({ dH, dS, entries: Object.freeze(entries) });
}

function calculateStemDuplexTm(
	stemThermodynamics,
	stemSequence,
	conditions,
	terminalMismatches = [],
) {
	const normalized = normalizeConditions(conditions);
	const terminal = sumTerminalMismatchParams(terminalMismatches);
	const salt = calculateOwczarzySaltCorrection({
		stemSequence,
		stemDeltaH: stemThermodynamics.dH,
		magnesiumMm: normalized.magnesiumMm,
		monovalentMm: normalized.monovalentMm,
	});
	const concentration = {
		concentrationUm: normalized.concentrationUm,
		limitingConcentrationUm: normalized.limitingConcentrationUm,
	};
	const withoutTerminalMismatches = calculateTmFromThermodynamics({
		dH: stemThermodynamics.dH,
		dS: stemThermodynamics.dS,
		saltCorrection: salt.saltCorrection,
		...concentration,
	});
	const withTerminalMismatches = calculateTmFromThermodynamics({
		dH: stemThermodynamics.dH + terminal.dH,
		dS: stemThermodynamics.dS + terminal.dS,
		saltCorrection: salt.saltCorrection,
		...concentration,
	});
	return Object.freeze({
		stemTm: withTerminalMismatches,
		stemTmWithoutTerminalMismatches: withoutTerminalMismatches,
		terminalMismatches: terminal,
		salt,
		conditions: normalized,
	});
}

/**
 * Wittwer/empirical snapback Tm.
 *
 * The empirical loop term always uses log10, matching the prior uSnapback web
 * app and the workbook compatibility calculation. An explicit non-log10
 * `conditions.wittwerLogBase` or legacy `options.logBase` value is rejected.
 * The optional fifth argument also permits the two real terminal mismatch
 * tetrads to be included in the stem Tm; callers using the legacy four-argument
 * signature calculate the core stem only.
 */
export function calculateSnapbackTmWittwer(
	stemSequence,
	loopLength,
	mismatch = undefined,
	conditions = {},
	options = {},
) {
	if (!Number.isInteger(loopLength) || loopLength < 3) {
		throw new Error('loopLength must be an integer of at least 3 bases.');
	}
	const normalizedConditions = normalizeConditions(conditions);
	const logBase = options.logBase ?? normalizedConditions.wittwerLogBase;
	const stem = calculateDuplexThermodynamics(stemSequence, mismatch);
	const duplex = calculateStemDuplexTm(
		stem,
		stemSequence,
		normalizedConditions,
		options.terminalMismatches ?? [],
	);
	const tm = calculateWittwerTmFromStemTm(
		duplex.stemTm,
		loopLength,
		logBase,
	);
	const result = Object.freeze({
		tm: roundTm(tm),
		unroundedTm: tm,
		loopLength,
		stem,
		duplex,
		logBase,
		coefficients: WITTWER_COEFFICIENTS,
	});
	return options.returnDetails ? result : result.tm;
}

export function getSantaLuciaHairpinLoopParams(loopLength) {
	if (!Number.isInteger(loopLength) || loopLength < 3) {
		throw new Error('SantaLucia/Hicks hairpin loops require an integer N >= 3.');
	}
	const exact = HAIRPIN_LOOP_ANCHORS[loopLength];
	if (exact) {
		return Object.freeze({
			...exact,
			N: loopLength,
			interpolated: false,
		});
	}

	let dG37;
	if (loopLength <= 30) {
		const sizes = Object.keys(HAIRPIN_LOOP_ANCHORS)
			.map(Number)
			.sort((a, b) => a - b);
		const upperIndex = sizes.findIndex((size) => size > loopLength);
		const lower = sizes[upperIndex - 1];
		const upper = sizes[upperIndex];
		const fraction = (loopLength - lower) / (upper - lower);
		dG37 =
			HAIRPIN_LOOP_ANCHORS[lower].dG37 +
			fraction *
				(HAIRPIN_LOOP_ANCHORS[upper].dG37 -
					HAIRPIN_LOOP_ANCHORS[lower].dG37);
	} else {
		const extrapolationCoefficient =
			(2.44 * GAS_CONSTANT_CAL * TEMPERATURE_37_K) / 1000;
		dG37 =
			HAIRPIN_LOOP_ANCHORS[30].dG37 +
			extrapolationCoefficient * Math.log(loopLength / 30);
	}

	return Object.freeze({
		N: loopLength,
		dG37,
		dH: 0,
		dS: (-dG37 * 1000) / TEMPERATURE_37_K,
		interpolated: loopLength <= 30,
	});
}

function validateExtendedSnapback(extended) {
	if (!extended || typeof extended !== 'object' || Array.isArray(extended)) {
		throw new Error('extendedSnapback must be an object.');
	}
	for (const key of [
		'fivePrimerLimSnapExtMismatches',
		'fivePrimeStem',
		'fivePrimeInnerLoopMismatches',
		'stuffBetween',
		'threePrimeInnerLoopMismatches',
		'threePrimeStem',
		'threePrimerLimSnapExtMismatches',
		'threePrimerRestOfAmplicon',
	]) {
		if (typeof extended[key] !== 'string' || !/^[ACGT]*$/.test(extended[key])) {
			throw new Error(`extendedSnapback.${key} must be an uppercase DNA string.`);
		}
	}
	validateDna('extendedSnapback.fivePrimeStem', extended.fivePrimeStem, 3);
	validateDna('extendedSnapback.threePrimeStem', extended.threePrimeStem, 3);
	if (extended.fivePrimeStem.length !== extended.threePrimeStem.length) {
		throw new Error('The two descriptive stem arms must have equal length.');
	}
	if (
		extended.fivePrimerLimSnapExtMismatches.length !== 1 ||
		extended.threePrimerLimSnapExtMismatches.length !== 1
	) {
		throw new Error(
			'The snapback must contain exactly one extension-blocking mismatch base on each strand.',
		);
	}
	if (
		extended.fivePrimeInnerLoopMismatches.length !==
			extended.threePrimeInnerLoopMismatches.length ||
		extended.fivePrimeInnerLoopMismatches.length > 1
	) {
		throw new Error(
			'The loop end must contain either one engineered mismatch pair or one natural mismatch already in the loop.',
		);
	}
	const snv = extended.snvOnThreePrimeStem;
	const tailSnv = extended.snvOnFivePrimeStem;
	if (
		!snv ||
		!Number.isInteger(snv.indexInThreePrimeStem) ||
		!VALID_DNA_BASES.has(snv.wildBase) ||
		!VALID_DNA_BASES.has(snv.variantBase) ||
		snv.wildBase === snv.variantBase ||
		!tailSnv ||
		!Number.isInteger(tailSnv.indexInFivePrimeStem) ||
		!VALID_DNA_BASES.has(tailSnv.tailBaseAtSNV)
	) {
		throw new Error('The extended snapback has invalid internal-SNV metadata.');
	}
	if (
		snv.indexInThreePrimeStem < MIN_INTERNAL_MATCHED_FLANK ||
		snv.indexInThreePrimeStem >=
			extended.threePrimeStem.length - MIN_INTERNAL_MATCHED_FLANK
	) {
		throw new Error(
			`The extended snapback SNV must have at least ${MIN_INTERNAL_MATCHED_FLANK} matched stem bases on each side.`,
		);
	}
	if (extended.threePrimeStem[snv.indexInThreePrimeStem] !== snv.wildBase) {
		throw new Error('The SNV wild base does not match threePrimeStem.');
	}
	const expectedTailIndex =
		extended.fivePrimeStem.length - 1 - snv.indexInThreePrimeStem;
	if (
		tailSnv.indexInFivePrimeStem !== expectedTailIndex ||
		extended.fivePrimeStem[expectedTailIndex] !== tailSnv.tailBaseAtSNV
	) {
		throw new Error('The tail SNV metadata does not match fivePrimeStem.');
	}
	const compWildBase = DNA_COMPLEMENT[snv.wildBase];
	const compVariantBase = DNA_COMPLEMENT[snv.variantBase];
	const matchesWild = tailSnv.tailBaseAtSNV === compWildBase;
	const matchesVariant = tailSnv.tailBaseAtSNV === compVariantBase;
	if (
		tailSnv.compWildBase !== compWildBase ||
		tailSnv.compVariantBase !== compVariantBase ||
		tailSnv.matchesWild !== matchesWild ||
		tailSnv.matchesVariant !== matchesVariant ||
		matchesWild === matchesVariant
	) {
		throw new Error('The tail allele-match metadata is internally inconsistent.');
	}
}

function getAlleleStructures(extended) {
	validateExtendedSnapback(extended);
	const snv = extended.snvOnThreePrimeStem;
	const index = snv.indexInThreePrimeStem;
	const wildTop = extended.threePrimeStem;
	const variantTop =
		wildTop.slice(0, index) + snv.variantBase + wildTop.slice(index + 1);
	const bottom = extended.fivePrimeStem;
	const wildStem = calculateDuplexThermodynamicsFromStrands(wildTop, bottom);
	const variantStem = calculateDuplexThermodynamicsFromStrands(
		variantTop,
		bottom,
	);
	if (
		wildStem.mismatchPositions.length +
			variantStem.mismatchPositions.length !==
		1
	) {
		throw new Error(
			'The fixed snapback tail must match exactly one of the wild or variant alleles.',
		);
	}
	return Object.freeze({
		wild: Object.freeze({ top: wildTop, stem: wildStem }),
		variant: Object.freeze({ top: variantTop, stem: variantStem }),
	});
}

function getStructureTerminalMismatches(extended) {
	const bottomAligned = reverse(extended.fivePrimeStem);
	const topStem = extended.threePrimeStem;

	const loopTopOutside =
		extended.threePrimeInnerLoopMismatches.slice(-1) ||
		extended.stuffBetween.slice(-1);
	const loopBottomOutside =
		extended.fivePrimeInnerLoopMismatches[0] || extended.stuffBetween[0];
	const tailTopOutside = extended.threePrimerLimSnapExtMismatches[0];
	const tailBottomOutside = extended.fivePrimerLimSnapExtMismatches.slice(-1);
	if (
		!loopTopOutside ||
		!loopBottomOutside ||
		!tailTopOutside ||
		!tailBottomOutside
	) {
		throw new Error(
			'Both the loop-side and extension-blocking terminal mismatch contexts are required.',
		);
	}
	if (isCanonicalPair(loopTopOutside, loopBottomOutside)) {
		throw new Error('The loop-side terminal pair is not a mismatch.');
	}
	if (isCanonicalPair(tailTopOutside, tailBottomOutside)) {
		throw new Error('The extension-blocking terminal pair is not a mismatch.');
	}

	return Object.freeze([
		Object.freeze({
			label: extended.fivePrimeInnerLoopMismatches
				? 'engineered-loop-side'
				: 'natural-loop-side',
			top2: `${loopTopOutside}${topStem[0]}`,
			bottom2: `${loopBottomOutside}${bottomAligned[0]}`,
		}),
		Object.freeze({
			label: 'extension-blocking-end',
			top2: `${topStem.slice(-1)}${tailTopOutside}`,
			bottom2: `${bottomAligned.slice(-1)}${tailBottomOutside}`,
		}),
	]);
}

function getStructureLoopLength(extended) {
	return (
		extended.stuffBetween.length +
		extended.threePrimeInnerLoopMismatches.length +
		extended.fivePrimeInnerLoopMismatches.length
	);
}

function calculateAlleleWittwer(
	allele,
	loopLength,
	terminalMismatches,
	conditions,
) {
	const duplex = calculateStemDuplexTm(
		allele.stem,
		allele.top,
		conditions,
		terminalMismatches,
	);
	const logBase = duplex.conditions.wittwerLogBase;
	const unroundedTm = calculateWittwerTmFromStemTm(
		duplex.stemTm,
		loopLength,
		logBase,
	);
	return Object.freeze({
		tm: roundTm(unroundedTm),
		unroundedTm,
		duplex,
		logBase,
	});
}

/** Wittwer/empirical calculation for a complete web-app snapback structure. */
export function calculateSnapbackTmWittwerFromStructure(
	extendedSnapback,
	conditions = {},
) {
	const alleles = getAlleleStructures(extendedSnapback);
	const terminalMismatches = getStructureTerminalMismatches(extendedSnapback);
	const loopLength = getStructureLoopLength(extendedSnapback);
	const wild = calculateAlleleWittwer(
		alleles.wild,
		loopLength,
		terminalMismatches,
		conditions,
	);
	const variant = calculateAlleleWittwer(
		alleles.variant,
		loopLength,
		terminalMismatches,
		conditions,
	);
	return Object.freeze({
		method: 'Wittwer/empirical',
		wildTm: wild.tm,
		variantTm: variant.tm,
		loopLength,
		logBase: wild.logBase,
		coefficients: WITTWER_COEFFICIENTS,
		terminalMismatches,
		alleles: Object.freeze({ wild, variant }),
	});
}

function calculateAlleleSantaLucia(
	allele,
	loop,
	terminalMismatches,
	conditions,
) {
	const normalized = normalizeConditions(conditions);
	const terminal = sumTerminalMismatchParams(terminalMismatches);
	const salt = calculateOwczarzySaltCorrection({
		stemSequence: allele.top,
		stemDeltaH: allele.stem.dH,
		magnesiumMm: normalized.magnesiumMm,
		monovalentMm: normalized.monovalentMm,
	});

	// The methods write-up defines a snapback by removing the ordinary duplex
	// initiation from the stem, then adding the loop and both terminal mismatches.
	const stemWithoutInitiation = Object.freeze({
		dH: allele.stem.dH - DUPLEX_INITIATION.dH,
		dS: allele.stem.dS - DUPLEX_INITIATION.dS,
	});
	const dH = stemWithoutInitiation.dH + loop.dH + terminal.dH;
	const dS = stemWithoutInitiation.dS + loop.dS + terminal.dS;
	const unroundedTm = calculateTmFromThermodynamics({
		dH,
		dS,
		saltCorrection: salt.saltCorrection,
	});

	return Object.freeze({
		tm: roundTm(unroundedTm),
		unroundedTm,
		dH,
		dS,
		stem: allele.stem,
		stemWithoutInitiation,
		terminalMismatches: terminal,
		salt,
	});
}

/**
 * Full SantaLucia/Hicks snapback Tm for the descriptive structure already
 * produced by the web app. This is intramolecular: no concentration term is
 * present in the final hairpin equation. Owczarzy is calculated from, and
 * applied to, the paired stem only.
 */
export function calculateSnapbackTmSantaLucia(
	extendedSnapback,
	conditions = {},
) {
	const alleles = getAlleleStructures(extendedSnapback);
	const terminalMismatches = getStructureTerminalMismatches(extendedSnapback);
	const loopLength = getStructureLoopLength(extendedSnapback);
	const loop = getSantaLuciaHairpinLoopParams(loopLength);
	const wild = calculateAlleleSantaLucia(
		alleles.wild,
		loop,
		terminalMismatches,
		conditions,
	);
	const variant = calculateAlleleSantaLucia(
		alleles.variant,
		loop,
		terminalMismatches,
		conditions,
	);

	return Object.freeze({
		method: 'SantaLucia/Hicks 2004 + Owczarzy 2008',
		wildTm: wild.tm,
		variantTm: variant.tm,
		components: Object.freeze({
			loop,
			terminalMismatches,
			stem: Object.freeze({
				wild: wild.stem,
				variant: variant.stem,
			}),
		}),
		sums: Object.freeze({
			wild: Object.freeze({
				dH: wild.dH,
				dS: wild.dS,
				saltCorrection: wild.salt.saltCorrection,
			}),
			variant: Object.freeze({
				dH: variant.dH,
				dS: variant.dS,
				saltCorrection: variant.salt.saltCorrection,
			}),
		}),
		alleles: Object.freeze({ wild, variant }),
	});
}

export { normalizeConditions as normalizeInHouseTmConditions };
