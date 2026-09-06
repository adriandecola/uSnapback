import { VALID_DNA_BASES } from './parameters.js';

const OWCZARZY_COEFFICIENTS = Object.freeze({
	a: 3.92e-5,
	b: -9.11e-6,
	c: 6.26e-5,
	d: 1.42e-5,
	e: -4.82e-4,
	f: 5.25e-4,
	g: 8.31e-5,
});

function validateSequence(sequence) {
	if (
		typeof sequence !== 'string' ||
		sequence.length < 2 ||
		![...sequence].every((base) => VALID_DNA_BASES.has(base))
	) {
		throw new Error('Stem sequence must contain at least two uppercase DNA bases.');
	}
}

function validateConcentration(name, value) {
	if (typeof value !== 'number' || !Number.isFinite(value) || value < 0) {
		throw new Error(`${name} must be a finite, non-negative number in mM.`);
	}
}

/**
 * Owczarzy et al. (2008) entropy-equivalent salt correction.
 *
 * The paper expresses the correction as a change in reciprocal Tm. Multiplying
 * that factor by the stem enthalpy in cal/mol yields the entropy-like term that
 * is added to a nearest-neighbour Tm denominator. Ionic inputs are accepted in
 * mM at the application boundary and converted to molar here.
 *
 * Only the paired stem contributes N, GC fraction, and dH. Loop and terminal-
 * mismatch thermodynamics are deliberately excluded from this salt term.
 */
export function calculateOwczarzySaltCorrection({
	stemSequence,
	stemDeltaH,
	magnesiumMm,
	monovalentMm,
}) {
	validateSequence(stemSequence);
	if (
		typeof stemDeltaH !== 'number' ||
		!Number.isFinite(stemDeltaH) ||
		stemDeltaH >= 0
	) {
		throw new Error('stemDeltaH must be a finite, negative value in kcal/mol.');
	}
	validateConcentration('magnesiumMm', magnesiumMm);
	validateConcentration('monovalentMm', monovalentMm);
	if (magnesiumMm === 0 && monovalentMm === 0) {
		throw new Error(
			'Owczarzy salt correction requires magnesiumMm or monovalentMm to be greater than zero.',
		);
	}

	const magnesiumM = magnesiumMm / 1000;
	const monovalentM = monovalentMm / 1000;
	const n = stemSequence.length;
	const gcFraction =
		[...stemSequence].filter((base) => base === 'G' || base === 'C').length /
		n;

	let regime;
	let ratio;
	let correctionFactor;
	let coefficients = null;

	if (magnesiumM === 0) {
		regime = 'monovalent-only';
		ratio = 0;
		const lnMonovalent = Math.log(monovalentM);
		correctionFactor =
			(4.29 * gcFraction - 3.95) * 1e-5 * lnMonovalent +
			9.4e-6 * lnMonovalent ** 2;
	} else if (monovalentM === 0) {
		regime = 'magnesium-only';
		ratio = Number.POSITIVE_INFINITY;
		coefficients = OWczarzyMagnesiumCoefficients(
			OWCZARZY_COEFFICIENTS,
			magnesiumM,
			gcFraction,
			n,
		);
		correctionFactor = coefficients.correctionFactor;
	} else {
		ratio = Math.sqrt(magnesiumM) / monovalentM;
		if (ratio < 0.22) {
			regime = 'monovalent-dominant';
			const lnMonovalent = Math.log(monovalentM);
			correctionFactor =
				(4.29 * gcFraction - 3.95) * 1e-5 * lnMonovalent +
				9.4e-6 * lnMonovalent ** 2;
		} else {
			const base = OWczarzyAdjustedCoefficients(monovalentM, ratio);
			regime = ratio < 6 ? 'mixed-ion' : 'magnesium-dominant';
			coefficients = OWczarzyMagnesiumCoefficients(
				base,
				magnesiumM,
				gcFraction,
				n,
			);
			correctionFactor = coefficients.correctionFactor;
		}
	}

	const saltCorrection = stemDeltaH * 1000 * correctionFactor;
	if (!Number.isFinite(saltCorrection)) {
		throw new Error('Owczarzy salt correction produced a non-finite result.');
	}

	return Object.freeze({
		saltCorrection,
		correctionFactor,
		regime,
		ratio,
		magnesiumM,
		monovalentM,
		gcFraction,
		stemLength: n,
		coefficients,
	});
}

function OWczarzyAdjustedCoefficients(monovalentM, ratio) {
	if (ratio >= 6) return OWCZARZY_COEFFICIENTS;

	const lnMonovalent = Math.log(monovalentM);
	return Object.freeze({
		...OWCZARZY_COEFFICIENTS,
		a:
			OWCZARZY_COEFFICIENTS.a *
			(0.843 - 0.352 * Math.sqrt(monovalentM) * lnMonovalent),
		d:
			OWCZARZY_COEFFICIENTS.d *
			(1.279 -
				4.03e-3 * lnMonovalent -
				8.03e-3 * lnMonovalent ** 2),
		g:
			OWCZARZY_COEFFICIENTS.g *
			(0.486 -
				0.258 * lnMonovalent +
				5.25e-3 * lnMonovalent ** 3),
	});
}

function OWczarzyMagnesiumCoefficients(base, magnesiumM, gcFraction, n) {
	const lnMagnesium = Math.log(magnesiumM);
	const correctionFactor =
		base.a +
		base.b * lnMagnesium +
		gcFraction * (base.c + base.d * lnMagnesium) +
		(base.e +
			base.f * lnMagnesium +
			base.g * lnMagnesium ** 2) /
			(2 * (n - 1));

	return Object.freeze({ ...base, correctionFactor });
}

export { OWczarzyAdjustedCoefficients, OWczarzyMagnesiumCoefficients };
