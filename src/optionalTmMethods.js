/*
File:           optionalTmMethods.js
Description:    Optional Santa Lucia and Rochester snapback Tm methods.
Author:         Adrian deCola
Relative Path:  uSnapback/src/optionalTmMethods.js
*/

/*****************************************************************************************/
/*************************************** Constants ***************************************/
/*****************************************************************************************/
const NUCLEOTIDE_COMPLEMENT = { A: 'T', T: 'A', C: 'G', G: 'C' };
const VALID_BASES = new Set(['A', 'T', 'C', 'G']);
const TM_DECIMAL_PLACES = 2;
const CONC = 0.5;
const LIMITING_CONC = 0.5;

//──────────────────────────────────────────────────────────────────────────//
// Rochester Hairpin Loop Initiation Parameters
// For loop size N > 30, we use the corrected large-N formula:
//
//   ΔH°(n) = -14.1
//   ΔS°(n) = -[ (14.1 + 6.5 + 0.1*(N - 30) ]*1000/310.15
//
//──────────────────────────────────────────────────────────────────────────//
const HAIRPIN_LOOP_PARAMETER_ROCHESTER = deepFreeze({
	3: { dH: -0.6, dS: -12.9 },
	4: { dH: -2.3, dS: -18.4 },
	5: { dH: -12.1, dS: -50.3 },
	6: { dH: -16.7, dS: -67.4 },
	7: { dH: -13.6, dS: -57.4 },
	8: { dH: -14.1, dS: -59.0 },
	9: { dH: -14.1, dS: -59.3 },
	10: { dH: -14.1, dS: -59.6 },
	11: { dH: -14.1, dS: -60.0 },
	12: { dH: -14.1, dS: -60.3 },
	13: { dH: -14.1, dS: -60.6 },
	14: { dH: -14.1, dS: -60.9 },
	15: { dH: -14.1, dS: -61.6 },
	16: { dH: -14.1, dS: -61.9 },
	17: { dH: -14.1, dS: -62.2 },
	18: { dH: -14.1, dS: -62.6 },
	19: { dH: -14.1, dS: -62.9 },
	20: { dH: -14.1, dS: -63.2 },
	21: { dH: -14.1, dS: -63.5 },
	22: { dH: -14.1, dS: -63.8 },
	23: { dH: -14.1, dS: -64.2 },
	24: { dH: -14.1, dS: -64.5 },
	25: { dH: -14.1, dS: -64.8 },
	26: { dH: -14.1, dS: -65.1 },
	27: { dH: -14.1, dS: -65.5 },
	28: { dH: -14.1, dS: -65.8 },
	29: { dH: -14.1, dS: -66.1 },
	30: { dH: -14.1, dS: -66.4 },
});

//──────────────────────────────────────────────────────────────────────────//
// SantaLucia & Hicks (2004) Hairpin Loop Initiation Parameters
// Paper assumes delta H = 0 for all loop initiations
// For loop size N ≤ 30 that not listed in the table, we do linear interpolation.
// For loop size N > 30, we use the corrected large-N formula:
//
//   ΔH°(n) = 0
//   ΔS°(n) = -[ 6.3 + 1.50*ln(n/30) ]*1000/310.15
//
//──────────────────────────────────────────────────────────────────────────//
const HAIRPIN_LOOP_PARAMETERS_SANTA_LUCIA_HICKS = deepFreeze({
	3: { dH: 0.0, dS: -11.3 },
	4: { dH: 0.0, dS: -11.3 },
	5: { dH: 0.0, dS: -10.6 },
	6: { dH: 0.0, dS: -12.9 },
	7: { dH: 0.0, dS: -13.5 },
	8: { dH: 0.0, dS: -14.3 },
	9: { dH: 0.0, dS: -14.5 },
	10: { dH: 0.0, dS: -14.6 },
	12: { dH: 0.0, dS: -16.1 },
	14: { dH: 0.0, dS: -16.4 },
	16: { dH: 0.0, dS: -17.1 },
	18: { dH: 0.0, dS: -17.7 },
	20: { dH: 0.0, dS: -18.4 },
	25: { dH: 0.0, dS: -19.7 },
	30: { dH: 0.0, dS: -20.3 },
});

/**
 * Bommarito (2000) dangling-end parameters.
 * Units: dH (kcal/mol), dS (cal/K/mol).
 *
 * Orientation rules:
 *   - 5′-dangling end of form XY: dangling base is X.
 *   - 3′-dangling end of form XY: dangling base is Y.
 */
const DANGLING_END_PARAMS = deepFreeze({
	fivePrime: {
		AA: { dH: 0.2, dS: 2.3 },
		AC: { dH: -6.3, dS: -17.1 },
		AG: { dH: -3.7, dS: -10.0 },
		AT: { dH: -2.9, dS: -7.6 },
		CA: { dH: 0.6, dS: 3.3 },
		CC: { dH: -4.4, dS: -12.6 },
		CG: { dH: -4.0, dS: -11.9 },
		CT: { dH: -4.1, dS: -13.0 },
		GA: { dH: -1.1, dS: -1.6 },
		GC: { dH: -5.1, dS: -14.0 },
		GG: { dH: -3.9, dS: -10.9 },
		GT: { dH: -4.2, dS: -15.0 },
		TA: { dH: -6.9, dS: -20.0 },
		TC: { dH: -4.0, dS: -10.9 },
		TG: { dH: -4.9, dS: -13.8 },
		TT: { dH: -0.2, dS: -0.5 },
	},
	threePrime: {
		TA: { dH: -0.7, dS: -0.8 },
		GA: { dH: -2.1, dS: -3.9 },
		CA: { dH: -5.9, dS: -16.5 },
		AA: { dH: -0.5, dS: -1.1 },
		TC: { dH: 4.4, dS: 7.3 },
		GC: { dH: -0.2, dS: -0.1 },
		CC: { dH: -2.6, dS: -7.4 },
		AC: { dH: 4.7, dS: 14.2 },
		TG: { dH: -1.6, dS: -3.6 },
		GG: { dH: -3.9, dS: -11.2 },
		CG: { dH: -3.2, dS: -10.4 },
		AG: { dH: -4.1, dS: -13.1 },
		TT: { dH: 2.9, dS: 10.4 },
		GT: { dH: -4.4, dS: -13.1 },
		CT: { dH: -5.2, dS: -15.0 },
		AT: { dH: -3.8, dS: -12.6 },
	},
});

/**
 * Orientation queries used by the lookup helpers.
 * Accepted inputs for 5′:  '5p', 'fivePrime'
 * Accepted inputs for 3′:  '3p', 'threePrime'
 */
const DANGLING_ORIENTATION = deepFreeze({
	FIVE_PRIME: 'fivePrime',
	THREE_PRIME: 'threePrime',
});

//──────────────────────────────────────────────────────────────────────────//
// Terminal Mismatches (Rochester)
// Units: dH (kcal/mol), dS (cal/K/mol)
// Key format: "TOP2/BOTTOM2"
//   • TOP2    is the 5'3' NN pair on the top strand
//   • BOTTOM2 is the 3'5' NN pair on the bottom strand
//
// Both printed forms from NNDB (“X/Y or Y'/X'”) are included as distinct keys
// mapping to the SAME numeric values.
//
// Remember to only use these parameters when the last base pair is a mismatch.
// Do not use this on penultimate mismatches.
//──────────────────────────────────────────────────────────────────────────//
const TERMINAL_MISMATCH_PARAMS = deepFreeze({
	// ── A/T neighbor match ────────────────────────────────────────────────
	'AA/TA': { dH: 4.0, dS: 15.2 },
	'AT/AA': { dH: 4.0, dS: 15.2 },

	'AA/TC': { dH: 5.2, dS: 17.7 },
	'CT/AA': { dH: 5.2, dS: 17.7 },

	'AA/TG': { dH: 4.7, dS: 16.8 },
	'GT/AA': { dH: 4.7, dS: 16.8 },

	'AC/TA': { dH: 4.0, dS: 14.8 },
	'AT/CA': { dH: 4.0, dS: 14.8 },

	'AC/TC': { dH: 5.2, dS: 17.4 },
	'CT/CA': { dH: 5.2, dS: 17.4 },

	'AC/TT': { dH: 4.3, dS: 14.8 },
	'TT/CA': { dH: 4.3, dS: 14.8 },

	'AG/TA': { dH: 4.0, dS: 14.8 },
	'AT/GA': { dH: 4.0, dS: 14.8 },

	'AG/TG': { dH: 4.7, dS: 16.4 },
	'GT/GA': { dH: 4.7, dS: 16.4 },

	'AG/TT': { dH: 4.3, dS: 15.5 },
	'TT/GA': { dH: 4.3, dS: 15.5 },

	'AT/TC': { dH: 5.2, dS: 17.7 },
	'CT/TA': { dH: 5.2, dS: 17.7 },

	'AT/TG': { dH: 4.7, dS: 16.8 },
	'GT/TA': { dH: 4.7, dS: 16.8 },

	'AT/TT': { dH: 4.3, dS: 15.2 },
	'TT/TA': { dH: 4.3, dS: 15.2 },

	// ── T/A neighbor match ────────────────────────────────────────────────
	'TA/AA': { dH: 2.3, dS: 9.4 },
	'AA/AT': { dH: 2.3, dS: 9.4 },

	'TA/AC': { dH: 1.7, dS: 6.8 },
	'CA/AT': { dH: 1.7, dS: 6.8 },

	'TA/AG': { dH: 0.7, dS: 3.9 },
	'GA/AT': { dH: 0.7, dS: 3.9 },

	'TC/AA': { dH: 2.3, dS: 9.0 },
	'AA/CT': { dH: 2.3, dS: 9.0 },

	'TC/AC': { dH: 1.7, dS: 6.1 },
	'CA/CT': { dH: 1.7, dS: 6.1 },

	'TC/AT': { dH: -5.8, dS: -17.1 },
	'TA/CT': { dH: -5.8, dS: -17.1 },

	'TG/AA': { dH: 2.3, dS: 9.4 },
	'AA/GT': { dH: 2.3, dS: 9.4 },

	'TG/AG': { dH: 0.7, dS: 3.5 },
	'GA/GT': { dH: 0.7, dS: 3.5 },

	'TG/AT': { dH: -5.8, dS: -17.1 },
	'TA/GT': { dH: -5.8, dS: -17.1 },

	'TT/AC': { dH: 1.7, dS: 6.4 },
	'CA/TT': { dH: 1.7, dS: 6.4 },

	'TT/AG': { dH: 0.7, dS: 4.2 },
	'GA/TT': { dH: 0.7, dS: 4.2 },

	'TT/AT': { dH: -5.8, dS: -17.7 },
	'TA/TT': { dH: -5.8, dS: -17.7 },

	// ── C/G neighbor match ────────────────────────────────────────────────
	'CA/GA': { dH: -4.6, dS: -11.6 },
	'AG/AC': { dH: -4.6, dS: -11.6 },

	'CA/GC': { dH: -4.2, dS: -11.0 },
	'CG/AC': { dH: -4.2, dS: -11.0 },

	'CA/GG': { dH: -4.7, dS: -12.3 },
	'GG/AC': { dH: -4.7, dS: -12.3 },

	'CC/GA': { dH: -4.6, dS: -12.3 },
	'AG/CC': { dH: -4.6, dS: -12.3 },

	'CC/GC': { dH: -4.2, dS: -11.9 },
	'CG/CC': { dH: -4.2, dS: -11.9 },

	'CC/GT': { dH: -5.0, dS: -13.9 },
	'TG/CC': { dH: -5.0, dS: -13.9 },

	'CG/GA': { dH: -4.6, dS: -11.6 },
	'AG/GC': { dH: -4.6, dS: -11.6 },

	'CG/GG': { dH: -4.7, dS: -12.3 },
	'GG/GC': { dH: -4.7, dS: -12.3 },

	'CG/GT': { dH: -5.0, dS: -12.9 },
	'TG/GC': { dH: -5.0, dS: -12.9 },

	'CT/GC': { dH: -4.2, dS: -11.6 },
	'CG/TC': { dH: -4.2, dS: -11.6 },

	'CT/GG': { dH: -4.7, dS: -12.3 },
	'GG/TC': { dH: -4.7, dS: -12.3 },

	'CT/GT': { dH: -5.0, dS: -13.2 },
	'TG/TC': { dH: -5.0, dS: -13.2 },

	// ── G/C neighbor match ────────────────────────────────────────────────
	'GA/CA': { dH: -6.2, dS: -16.8 },
	'AC/AG': { dH: -6.2, dS: -16.8 },

	'GA/CC': { dH: -3.3, dS: -8.4 },
	'CC/AG': { dH: -3.3, dS: -8.4 },

	'GA/CG': { dH: -5.1, dS: -13.9 },
	'GC/AG': { dH: -5.1, dS: -13.9 },

	'GC/CA': { dH: -6.2, dS: -16.8 },
	'AC/CG': { dH: -6.2, dS: -16.8 },

	'GC/CC': { dH: -3.3, dS: -8.7 },
	'CC/CG': { dH: -3.3, dS: -8.7 },

	'GC/CT': { dH: -0.9, dS: -0.6 },
	'TC/CG': { dH: -0.9, dS: -0.6 },

	'GG/CA': { dH: -6.2, dS: -16.8 },
	'AC/GG': { dH: -6.2, dS: -16.8 },

	'GG/CG': { dH: -5.1, dS: -13.2 },
	'GC/GG': { dH: -5.1, dS: -13.2 },

	'GG/CT': { dH: -0.9, dS: -0.3 },
	'TC/GG': { dH: -0.9, dS: -0.3 },

	'GT/CC': { dH: -3.3, dS: -8.7 },
	'CC/TG': { dH: -3.3, dS: -8.7 },

	'GT/CG': { dH: -5.1, dS: -13.5 },
	'GC/TG': { dH: -5.1, dS: -13.5 },

	'GT/CT': { dH: -0.9, dS: 0.0 },
	'TC/TG': { dH: -0.9, dS: 0.0 },
});


/*****************************************************************************************/
/************************ Optional Santa Lucia / Rochester Methods ************************/
/*****************************************************************************************/
/**
 * Calculates snapback Tm (wild vs variant) for the extended snapback using:
 *   1) Rochester hairpin-loop initiation parameters (loop size = stuffBetween + both strands of the optional inner-loop mismatch)
 *   2) A 5′-dangling-end correction on the loop side of the stem
 *   3) Designed terminal-mismatch correction on the 3′ end of the stem
 *   4) Stem nearest-neighbor thermodynamics with/without the SNV mismatch
 *
 * The “wild” vs “variant” difference is applied only at the stem via getThermoParams:
 *   - If descriptiveExtendedSnapback.snvOnFivePrimeStem.matchesWild === true,
 *     then the WILD case uses the matched stem (no mismatch object),
 *     and the VARIANT case uses the mismatch (position=indexInThreePrimeStem, type=tailBaseAtSNV).
 *   - If matchesVariant === true, the roles are reversed.
 *
 * Concentrations are in µM to match getThermoParams; they default to global constants if provided.
 *
 * Returns per-allele Tm in °C plus a detailed breakdown of the summed ΔH/ΔS components.
 *
 * Assumptions
 * - All strings on descriptiveExtendedSnapback are 5'→3' and already validated upstream.\
 * - Designed terminal-mismatch correction is applied once, using the last paired base on threePrimeStem and the
 *   first mismatch base on each strand (threePrimerLimSnapExtMismatches[0] and fivePrimerLimSnapExtMismatches[0]).
 *
 * Parameters
 * @param {DescriptiveExtendedSnapback}  extended
 *
 * @returns {Promise<{
 *   wildTm: number,
 *   variantTm: number,
 *   components: {
 *     loop: { N: number, dH: number, dS: number, model: 'Rochester' },
 *     dangling5p: null | { step: string, dH: number, dS: number },
 *     terminalMismatch3p: null | { top2: string, bottom2: string, dH: number, dS: number },
 *     stem: {
 *       matched: { dH: number, dS: number, saltCorrection: number, mismatch: null },
 *       mismatched: { dH: number, dS: number, saltCorrection: number, mismatch: { position: number, type: string } }
 *     }
 *   },
 *   sums: {
 *     wild: { dH: number, dS: number, saltCorrection: number },
 *     variant: { dH: number, dS: number, saltCorrection: number }
 *   }
 * }>}
 *
 * Throws
 * - If extended object is malformed or any sequence constraint would make indices out of bounds.
 */
async function calculateSnapbackTmRochester(extended, options = {}) {
	//──────────────────────────────────────────────────────────────────────────//
	//                            Parameter Checking                            //
	//──────────────────────────────────────────────────────────────────────────//

	if (typeof extended !== 'object' || !extended) {
		throw new Error('extended must be a non-null object.');
	}

	const {
		getThermoParams,
		conc = CONC,
		limitingConc = LIMITING_CONC,
	} = options;
	if (typeof getThermoParams !== 'function') {
		throw new Error('calculateSnapbackTmRochester requires options.getThermoParams.');
	}

	const requiredStrings = [
		'fivePrimerLimSnapExtMismatches',
		'fivePrimeStem',
		'fivePrimeInnerLoopMismatches',
		'stuffBetween',
		'threePrimeInnerLoopMismatches',
		'threePrimeStem',
		'threePrimerLimSnapExtMismatches',
		'threePrimerRestOfAmplicon',
	];
	for (const key of requiredStrings) {
		if (
			typeof extended[key] !== 'string' ||
			!/^[ACGT]*$/.test(extended[key])
		) {
			throw new Error(
				`extended.${key} must be an uppercase DNA string (A/T/C/G).`,
			);
		}
	}
	if (
		typeof extended.snvOnThreePrimeStem !== 'object' ||
		extended.snvOnThreePrimeStem === null ||
		!Number.isInteger(extended.snvOnThreePrimeStem.indexInThreePrimeStem) ||
		extended.snvOnThreePrimeStem.indexInThreePrimeStem < 0 ||
		typeof extended.snvOnThreePrimeStem.wildBase !== 'string' ||
		typeof extended.snvOnThreePrimeStem.variantBase !== 'string'
	) {
		throw new Error(
			'extended.snvOnThreePrimeStem must be { indexInThreePrimeStem:int≥0, wildBase:str, variantBase:str }.',
		);
	}
	if (
		typeof extended.snvOnFivePrimeStem !== 'object' ||
		extended.snvOnFivePrimeStem === null ||
		!Number.isInteger(extended.snvOnFivePrimeStem.indexInFivePrimeStem) ||
		extended.snvOnFivePrimeStem.indexInFivePrimeStem < 0 ||
		typeof extended.snvOnFivePrimeStem.tailBaseAtSNV !== 'string' ||
		typeof extended.snvOnFivePrimeStem.compWildBase !== 'string' ||
		typeof extended.snvOnFivePrimeStem.compVariantBase !== 'string'
	) {
		throw new Error(
			'extended.snvOnFivePrimeStem must include indexInFivePrimeStem:int≥0, tailBaseAtSNV, compWildBase, compVariantBase.',
		);
	}

	// Aliases
	const fivePrimerLimSnapExtMismatches =
		extended.fivePrimerLimSnapExtMismatches;
	const fivePrimeStem = extended.fivePrimeStem;
	const fivePrimeInnerLoopMismatches = extended.fivePrimeInnerLoopMismatches;
	const stuffBetween = extended.stuffBetween;
	const threePrimeInnerLoopMismatches =
		extended.threePrimeInnerLoopMismatches;
	const threePrimeStem = extended.threePrimeStem;
	const threePrimerLimSnapExtMismatches =
		extended.threePrimerLimSnapExtMismatches;

	const snvIdx = extended.snvOnThreePrimeStem.indexInThreePrimeStem;
	const wildBaseAtSNV = extended.snvOnThreePrimeStem.wildBase;
	const variantBaseAtSNV = extended.snvOnThreePrimeStem.variantBase;
	const tailBaseAtSNV = extended.snvOnFivePrimeStem.tailBaseAtSNV;
	const matchesWild = extended.snvOnFivePrimeStem.matchesWild;
	const matchesVariant = extended.snvOnFivePrimeStem.matchesVariant;

	if (snvIdx >= threePrimeStem.length) {
		throw new Error(
			`SNV index ${snvIdx} is out of bounds for threePrimeStem length ${threePrimeStem.length}.`,
		);
	}
	if (!VALID_BASES.has(wildBaseAtSNV) || !VALID_BASES.has(variantBaseAtSNV)) {
		throw new Error(
			'snvOnThreePrimeStem.wildBase and .variantBase must each be one of A/T/C/G.',
		);
	}
	if (!VALID_BASES.has(tailBaseAtSNV)) {
		throw new Error(
			'snvOnFivePrimeStem.tailBaseAtSNV must be one of A/T/C/G.',
		);
	}
	if (matchesWild === matchesVariant) {
		throw new Error(
			'Exactly one of snvOnFivePrimeStem.matchesWild and matchesVariant must be true.',
		);
	}
	if (threePrimeStem[snvIdx] !== wildBaseAtSNV) {
		throw new Error(
			`SNV metadata mismatch: threePrimeStem[${snvIdx}]="${threePrimeStem[snvIdx]}" does not match snvOnThreePrimeStem.wildBase="${wildBaseAtSNV}".`,
		);
	}

	const threePrimeStemVariant =
		threePrimeStem.slice(0, snvIdx) +
		variantBaseAtSNV +
		threePrimeStem.slice(snvIdx + 1);

	//──────────────────────────────────────────────────────────────────────────//
	// 1) Loop initiation (Rochester): N = stuffBetween + both strands of the optional inner-loop mismatch
	//──────────────────────────────────────────────────────────────────────────//
	const loopN = stuffBetween.length + 2 * fivePrimeInnerLoopMismatches.length;
	const loopParams = getRochesterHairpinLoopParams(loopN); // { dH (kcal/mol), dS (cal/K/mol) }

	// 2) Terminal mismatch at the 5′ end of the stem (loop side)
	//    Use the base immediately outside the stem (loop) and the first base inside the stem.
	let terminal5p = null;
	if (
		threePrimeStem.length >= 1 &&
		threePrimeInnerLoopMismatches.length >= 1 &&
		fivePrimeInnerLoopMismatches.length >= 1
	) {
		const topOutsideLeft =
			threePrimeInnerLoopMismatches[
				threePrimeInnerLoopMismatches.length - 1
			]; // from seq, 5' of stem
		const topFirstInStem = threePrimeStem[0];

		const bottomOutsideLeft = fivePrimeInnerLoopMismatches[0]; // 3' end of TM
		const bottomFirstInStem = NUCLEOTIDE_COMPLEMENT[topFirstInStem];

		const top2_left = `${topOutsideLeft}${topFirstInStem}`; // 5'→3' on top
		const bottom2_left = `${bottomOutsideLeft}${bottomFirstInStem}`; // 3'→5' on bottom

		const tm5 = getTerminalMismatchParams(top2_left, bottom2_left); // { dH, dS }
		terminal5p = {
			top2: top2_left,
			bottom2: bottom2_left,
			dH: tm5.dH,
			dS: tm5.dS,
		};
	}

	// 3) Terminal mismatch at the 3′ end of the stem (right side)
	//    Use the last base inside the stem and the base immediately outside (right).
	let terminal3p = null;
	if (
		threePrimeStem.length >= 1 &&
		threePrimerLimSnapExtMismatches.length >= 1 &&
		fivePrimerLimSnapExtMismatches.length >= 1
	) {
		const topLastInStem = threePrimeStem[threePrimeStem.length - 1];
		const topOutsideRight = threePrimerLimSnapExtMismatches[0];

		const bottomLastInStem = NUCLEOTIDE_COMPLEMENT[topLastInStem];
		const bottomOutsideRight =
			fivePrimerLimSnapExtMismatches[
				fivePrimerLimSnapExtMismatches.length - 1
			];

		const top2_right = `${topLastInStem}${topOutsideRight}`; // 5'→3' on top
		const bottom2_right = `${bottomLastInStem}${bottomOutsideRight}`; // 3'→5' on bottom

		const tm3 = getTerminalMismatchParams(top2_right, bottom2_right); // { dH, dS }
		terminal3p = {
			top2: top2_right,
			bottom2: bottom2_right,
			dH: tm3.dH,
			dS: tm3.dS,
		};
	}

	// 4) Compute allele-specific stem thermodynamics.
	//    The wild/variant difference comes from:
	//      - the top-strand base at SNV (wild stem vs variant stem sequence)
	//      - whether the fixed tail base at SNV is a match or mismatch for that allele
	const wildMismatchSpec = matchesWild
		? undefined
		: { position: snvIdx, type: tailBaseAtSNV };
	const variantMismatchSpec = matchesWild
		? { position: snvIdx, type: tailBaseAtSNV }
		: undefined;

	const [wildStem, variantStem] = await Promise.all([
		getThermoParams(threePrimeStem, conc, limitingConc, wildMismatchSpec),
		getThermoParams(
			threePrimeStemVariant,
			conc,
			limitingConc,
			variantMismatchSpec,
		),
	]);

	// Backward-compatible breakdown buckets.
	const stemMatched = matchesWild ? wildStem : variantStem;
	const stemMismatched = matchesWild ? variantStem : wildStem;

	const common_dH =
		loopParams.dH +
		(terminal5p ? terminal5p.dH : 0) +
		(terminal3p ? terminal3p.dH : 0);
	const common_dS =
		loopParams.dS +
		(terminal5p ? terminal5p.dS : 0) +
		(terminal3p ? terminal3p.dS : 0);

	const wildSum_dH = common_dH + wildStem.dH;
	const wildSum_dS = common_dS + wildStem.dS;
	const wildSalt = wildStem.saltCorrection || 0;

	const variantSum_dH = common_dH + variantStem.dH;
	const variantSum_dS = common_dS + variantStem.dS;
	const variantSalt = variantStem.saltCorrection || 0;

	// 6) Compute Tm (°C)
	const wildTm = calculateTm(
		wildSum_dH,
		wildSum_dS,
		undefined,
		undefined,
		false,
		wildSalt,
	);
	const variantTm = calculateTm(
		variantSum_dH,
		variantSum_dS,
		undefined,
		undefined,
		false,
		variantSalt,
	);

	// 7) Return breakdown and totals
	return {
		wildTm,
		variantTm,
		components: {
			loop: {
				N: loopN,
				dH: loopParams.dH,
				dS: loopParams.dS,
				model: 'Rochester',
			},
			terminalMismatch5p: terminal5p
				? {
						top2: terminal5p.top2,
						bottom2: terminal5p.bottom2,
						dH: terminal5p.dH,
						dS: terminal5p.dS,
					}
				: null,
			terminalMismatch3p: terminal3p
				? {
						top2: terminal3p.top2,
						bottom2: terminal3p.bottom2,
						dH: terminal3p.dH,
						dS: terminal3p.dS,
					}
				: null,
			stem: {
				matched: {
					dH: stemMatched.dH,
					dS: stemMatched.dS,
					saltCorrection: stemMatched.saltCorrection || 0,
					mismatch: null,
				},
				mismatched: {
					dH: stemMismatched.dH,
					dS: stemMismatched.dS,
					saltCorrection: stemMismatched.saltCorrection || 0,
					mismatch: {
						position: snvIdx,
						type: tailBaseAtSNV,
					},
				},
			},
		},
		sums: {
			wild: { dH: wildSum_dH, dS: wildSum_dS, saltCorrection: wildSalt },
			variant: {
				dH: variantSum_dH,
				dS: variantSum_dS,
				saltCorrection: variantSalt,
			},
		},
	};
}

/**
 * Calculates snapback Tm (wild vs variant) for the extended snapback using:
 *   1) SantaLucia hairpin-loop initiation parameters (loop size = stuffBetween + both strands of the optional inner-loop mismatch)
 *   2) A 5′-dangling-end correction on the loop side of the stem
 *   3) Designed terminal-mismatch correction on the 3′ end of the stem
 *   4) Stem nearest-neighbor thermodynamics with/without the SNV mismatch
 *
 * The “wild” vs “variant” difference is applied only at the stem via getThermoParams:
 *   - If descriptiveExtendedSnapback.snvOnFivePrimeStem.matchesWild === true,
 *     then the WILD case uses the matched stem (no mismatch object),
 *     and the VARIANT case uses the mismatch (position=indexInThreePrimeStem, type=tailBaseAtSNV).
 *   - If matchesVariant === true, the roles are reversed.
 *
 * Concentrations are in µM to match getThermoParams; they default to global constants if provided.
 *
 * Returns per-allele Tm in °C plus a detailed breakdown of the summed ΔH/ΔS components.
 *
 * Assumptions
 * - All strings on descriptiveExtendedSnapback are 5'→3' and already validated upstream.\
 * - Designed terminal-mismatch correction is applied once, using the last paired base on threePrimeStem and the
 *   first mismatch base on each strand (threePrimerLimSnapExtMismatches[0] and fivePrimerLimSnapExtMismatches[0]).
 *
 * Parameters
 * @param {DescriptiveExtendedSnapback}  extended
 * @property {string} fivePrimerLimSnapExtMismatches  Strong mismatch placed to prevent extension on the complementary snapback.
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
 * @returns {Promise<{
 *   wildTm: number,
 *   variantTm: number,
 *   components: {
 *     loop: { N: number, dH: number, dS: number, model: 'SantaLuciaHicks' },
 *     dangling5p: null | { step: string, dH: number, dS: number },
 *     terminalMismatch3p: null | { top2: string, bottom2: string, dH: number, dS: number },
 *     stem: {
 *       matched: { dH: number, dS: number, saltCorrection: number, mismatch: null },
 *       mismatched: { dH: number, dS: number, saltCorrection: number, mismatch: { position: number, type: string } }
 *     }
 *   },
 *   sums: {
 *     wild: { dH: number, dS: number, saltCorrection: number },
 *     variant: { dH: number, dS: number, saltCorrection: number }
 *   }
 * }>}
 *
 * Throws
 * - If extended object is malformed or any sequence constraint would make indices out of bounds.
 */
async function calculateSnapbackTmSantaLucia(extended, options = {}) {
	//──────────────────────────────────────────────────────────────────────────//
	//                            Parameter Checking                            //
	//──────────────────────────────────────────────────────────────────────────//
	if (typeof extended !== 'object' || !extended) {
		throw new Error('extended must be a non-null object.');
	}

	const {
		getThermoParams,
		conc = CONC,
		limitingConc = LIMITING_CONC,
	} = options;
	if (typeof getThermoParams !== 'function') {
		throw new Error('calculateSnapbackTmSantaLucia requires options.getThermoParams.');
	}

	const requiredStrings = [
		'fivePrimerLimSnapExtMismatches',
		'fivePrimeStem',
		'fivePrimeInnerLoopMismatches',
		'stuffBetween',
		'threePrimeInnerLoopMismatches',
		'threePrimeStem',
		'threePrimerLimSnapExtMismatches',
		'threePrimerRestOfAmplicon',
	];
	for (const key of requiredStrings) {
		if (
			typeof extended[key] !== 'string' ||
			!/^[ACGT]*$/.test(extended[key])
		) {
			throw new Error(
				`extended.${key} must be an uppercase DNA string (A/T/C/G).`,
			);
		}
	}
	if (
		typeof extended.snvOnThreePrimeStem !== 'object' ||
		extended.snvOnThreePrimeStem === null ||
		!Number.isInteger(extended.snvOnThreePrimeStem.indexInThreePrimeStem) ||
		extended.snvOnThreePrimeStem.indexInThreePrimeStem < 0 ||
		typeof extended.snvOnThreePrimeStem.wildBase !== 'string' ||
		typeof extended.snvOnThreePrimeStem.variantBase !== 'string'
	) {
		throw new Error(
			'extended.snvOnThreePrimeStem must be { indexInThreePrimeStem:int≥0, wildBase:str, variantBase:str }.',
		);
	}
	if (
		typeof extended.snvOnFivePrimeStem !== 'object' ||
		extended.snvOnFivePrimeStem === null ||
		!Number.isInteger(extended.snvOnFivePrimeStem.indexInFivePrimeStem) ||
		extended.snvOnFivePrimeStem.indexInFivePrimeStem < 0 ||
		typeof extended.snvOnFivePrimeStem.tailBaseAtSNV !== 'string' ||
		typeof extended.snvOnFivePrimeStem.compWildBase !== 'string' ||
		typeof extended.snvOnFivePrimeStem.compVariantBase !== 'string'
	) {
		throw new Error(
			'extended.snvOnFivePrimeStem must include indexInFivePrimeStem:int≥0, tailBaseAtSNV, compWildBase, compVariantBase.',
		);
	}

	// Aliases
	const fivePrimerLimSnapExtMismatches =
		extended.fivePrimerLimSnapExtMismatches;
	const fivePrimeStem = extended.fivePrimeStem;
	const fivePrimeInnerLoopMismatches = extended.fivePrimeInnerLoopMismatches;
	const stuffBetween = extended.stuffBetween;
	const threePrimeInnerLoopMismatches =
		extended.threePrimeInnerLoopMismatches;
	const threePrimeStem = extended.threePrimeStem;
	const threePrimerLimSnapExtMismatches =
		extended.threePrimerLimSnapExtMismatches;

	const snvIdx = extended.snvOnThreePrimeStem.indexInThreePrimeStem;
	const wildBaseAtSNV = extended.snvOnThreePrimeStem.wildBase;
	const variantBaseAtSNV = extended.snvOnThreePrimeStem.variantBase;
	const tailBaseAtSNV = extended.snvOnFivePrimeStem.tailBaseAtSNV;
	const matchesWild = extended.snvOnFivePrimeStem.matchesWild;
	const matchesVariant = extended.snvOnFivePrimeStem.matchesVariant;

	if (snvIdx >= threePrimeStem.length) {
		throw new Error(
			`SNV index ${snvIdx} is out of bounds for threePrimeStem length ${threePrimeStem.length}.`,
		);
	}
	if (!VALID_BASES.has(wildBaseAtSNV) || !VALID_BASES.has(variantBaseAtSNV)) {
		throw new Error(
			'snvOnThreePrimeStem.wildBase and .variantBase must each be one of A/T/C/G.',
		);
	}
	if (!VALID_BASES.has(tailBaseAtSNV)) {
		throw new Error(
			'snvOnFivePrimeStem.tailBaseAtSNV must be one of A/T/C/G.',
		);
	}
	if (matchesWild === matchesVariant) {
		throw new Error(
			'Exactly one of snvOnFivePrimeStem.matchesWild and matchesVariant must be true.',
		);
	}
	if (threePrimeStem[snvIdx] !== wildBaseAtSNV) {
		throw new Error(
			`SNV metadata mismatch: threePrimeStem[${snvIdx}]="${threePrimeStem[snvIdx]}" does not match snvOnThreePrimeStem.wildBase="${wildBaseAtSNV}".`,
		);
	}

	const threePrimeStemVariant =
		threePrimeStem.slice(0, snvIdx) +
		variantBaseAtSNV +
		threePrimeStem.slice(snvIdx + 1);

	//──────────────────────────────────────────────────────────────────────────//
	// 1) Loop initiation (SantaLucia-Hicks): N = stuffBetween + both strands of the optional inner-loop mismatch
	//──────────────────────────────────────────────────────────────────────────//
	const loopN = stuffBetween.length + 2 * fivePrimeInnerLoopMismatches.length;
	const loopParams = getSantaLuciaHicksHairpinParams(loopN); // { dH (kcal/mol), dS (cal/K/mol) }

	// 2) Terminal mismatch at the 5′ end of the stem (loop side)
	//    Use the base immediately outside the stem (loop) and the first base inside the stem.
	let terminal5p = null;
	if (
		threePrimeStem.length >= 1 &&
		threePrimeInnerLoopMismatches.length >= 1 &&
		fivePrimeInnerLoopMismatches.length >= 1
	) {
		const topOutsideLeft =
			threePrimeInnerLoopMismatches[
				threePrimeInnerLoopMismatches.length - 1
			]; // from seq, 5' of stem
		const topFirstInStem = threePrimeStem[0];

		const bottomOutsideLeft = fivePrimeInnerLoopMismatches[0]; // 3' end of TM
		const bottomFirstInStem = NUCLEOTIDE_COMPLEMENT[topFirstInStem];

		const top2_left = `${topOutsideLeft}${topFirstInStem}`; // 5'→3' on top
		const bottom2_left = `${bottomOutsideLeft}${bottomFirstInStem}`; // 3'→5' on bottom

		const tm5 = getTerminalMismatchParams(top2_left, bottom2_left); // { dH, dS }
		terminal5p = {
			top2: top2_left,
			bottom2: bottom2_left,
			dH: tm5.dH,
			dS: tm5.dS,
		};
	}

	// 3) Terminal mismatch at the 3′ end of the stem (right side)
	//    Use the last base inside the stem and the base immediately outside (right).
	let terminal3p = null;
	if (
		threePrimeStem.length >= 1 &&
		threePrimerLimSnapExtMismatches.length >= 1 &&
		fivePrimerLimSnapExtMismatches.length >= 1
	) {
		const topLastInStem = threePrimeStem[threePrimeStem.length - 1];
		const topOutsideRight = threePrimerLimSnapExtMismatches[0];

		const bottomLastInStem = NUCLEOTIDE_COMPLEMENT[topLastInStem];
		const bottomOutsideRight =
			fivePrimerLimSnapExtMismatches[
				fivePrimerLimSnapExtMismatches.length - 1
			];

		const top2_right = `${topLastInStem}${topOutsideRight}`; // 5'→3' on top
		const bottom2_right = `${bottomLastInStem}${bottomOutsideRight}`; // 3'→5' on bottom

		const tm3 = getTerminalMismatchParams(top2_right, bottom2_right); // { dH, dS }
		terminal3p = {
			top2: top2_right,
			bottom2: bottom2_right,
			dH: tm3.dH,
			dS: tm3.dS,
		};
	}

	// 4) Compute allele-specific stem thermodynamics.
	//    The wild/variant difference comes from:
	//      - the top-strand base at SNV (wild stem vs variant stem sequence)
	//      - whether the fixed tail base at SNV is a match or mismatch for that allele
	const wildMismatchSpec = matchesWild
		? undefined
		: { position: snvIdx, type: tailBaseAtSNV };
	const variantMismatchSpec = matchesWild
		? { position: snvIdx, type: tailBaseAtSNV }
		: undefined;

	const [wildStem, variantStem] = await Promise.all([
		getThermoParams(threePrimeStem, conc, limitingConc, wildMismatchSpec),
		getThermoParams(
			threePrimeStemVariant,
			conc,
			limitingConc,
			variantMismatchSpec,
		),
	]);

	// Backward-compatible breakdown buckets.
	const stemMatched = matchesWild ? wildStem : variantStem;
	const stemMismatched = matchesWild ? variantStem : wildStem;

	const common_dH =
		loopParams.dH +
		(terminal5p ? terminal5p.dH : 0) +
		(terminal3p ? terminal3p.dH : 0);
	const common_dS =
		loopParams.dS +
		(terminal5p ? terminal5p.dS : 0) +
		(terminal3p ? terminal3p.dS : 0);

	const wildSum_dH = common_dH + wildStem.dH;
	const wildSum_dS = common_dS + wildStem.dS;
	const wildSalt = wildStem.saltCorrection || 0;

	const variantSum_dH = common_dH + variantStem.dH;
	const variantSum_dS = common_dS + variantStem.dS;
	const variantSalt = variantStem.saltCorrection || 0;

	// 6) Compute Tm (°C)
	const wildTm = calculateTm(
		wildSum_dH,
		wildSum_dS,
		undefined,
		undefined,
		false,
		wildSalt,
	);
	const variantTm = calculateTm(
		variantSum_dH,
		variantSum_dS,
		undefined,
		undefined,
		false,
		variantSalt,
	);

	// 7) Return breakdown and totals
	return {
		wildTm,
		variantTm,
		components: {
			loop: {
				N: loopN,
				dH: loopParams.dH,
				dS: loopParams.dS,
				model: 'SantaLuciaHicks',
			},
			terminalMismatch5p: terminal5p
				? {
						top2: terminal5p.top2,
						bottom2: terminal5p.bottom2,
						dH: terminal5p.dH,
						dS: terminal5p.dS,
					}
				: null,
			terminalMismatch3p: terminal3p
				? {
						top2: terminal3p.top2,
						bottom2: terminal3p.bottom2,
						dH: terminal3p.dH,
						dS: terminal3p.dS,
					}
				: null,
			stem: {
				matched: {
					dH: stemMatched.dH,
					dS: stemMatched.dS,
					saltCorrection: stemMatched.saltCorrection || 0,
					mismatch: null,
				},
				mismatched: {
					dH: stemMismatched.dH,
					dS: stemMismatched.dS,
					saltCorrection: stemMismatched.saltCorrection || 0,
					mismatch: {
						position: snvIdx,
						type: tailBaseAtSNV,
					},
				},
			},
		},
		sums: {
			wild: { dH: wildSum_dH, dS: wildSum_dS, saltCorrection: wildSalt },
			variant: {
				dH: variantSum_dH,
				dS: variantSum_dS,
				saltCorrection: variantSalt,
			},
		},
	};
}

/**
 * Calculates the duplex melting temperature (Tm) at which 50% of the nucleic acid
 * is in the double-stranded state, using standard thermodynamic relations:
 *
 *      Tm = ΔH° / (ΔS° + saltCorrection - R·ln([AB]/([A][B])))
 *         = ΔH° / (ΔS° + saltCorrection + R·ln(([A][B])/[AB]))
 *
 * Concentration term at Tm (50% duplex):
 *  - General case (A ≠ B): let [AB]_50 = min([A],[B]) / 2.
 *        [A]_50 = [A] - [AB]_50
 *        [B]_50 = [B] - [AB]_50
 *        concentrationTerm = ([A]_50 · [B]_50) / [AB]_50
 *
 *    This reduces to the well-known special cases:
 *      • Equimolar, [A] = [B] = C0:
 *            [AB]_50 = C0/2,  [A]_50 = [B]_50 = C0/2
 *            ([A]_50 · [B]_50)/[AB]_50 = (C0/2 · C0/2)/(C0/2) = C0/2 = (Ct / 4)
 *            where Ct = [A] + [B] = 2·C0
 * 			  K = [A]/2
 *
 *      • Limiting strand B < A:
 *            [AB]_50 = [B]/2 ⇒ K = [A] - [B]/2
 * 					(this case works even is [A]=[B])
 *
 *  - Self-complementary (A + A ↔ AA):
 *        At Tm: [A]_50 = Ct/2, [AA]_50 = Ct/4  ⇒  ([A]_50·[A]_50)/[AA]_50 = Ct
 *        We expose this via a boolean flag `selfComplementary = true` and use:
 *            concentrationTerm = [A]
 *
 *
 *
 * Units and constants:
 *  - ΔH° must be provided in kcal/mol (internally converted to cal/mol).
 *  - ΔS° must be provided in cal/(K·mol).
 *  - saltCorrection is an additive entropy-like term in cal/(K·mol).
 *  - R (gas constant) = 1.98720425864 cal/(K·mol).
 *  - Concentrations are in micro Moles
 *
 * Returns:
 *  - Tm in degrees Celsius, rounded to TM_DECIMAL_PLACES.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Assumptions
 * ──────────────────────────────────────────────────────────────────────────
 * - sumDeltaH is the total enthalpy change (kcal/mol) for the duplex of interest,
 *   e.g., summed nearest-neighbor + correction terms (dangling ends, terminal mismatches, etc.).
 * - sumDeltaS is the total entropy change (cal/K/mol) for the same duplex.
 * - concA and concB are initial single-strand concentrations (uM).
 * - For self-complementary duplexes, pass selfComplementary = true and provide concA (> 0).
 *   concB is ignored in that case. Only do this if the DNA is self-complimentary AND there
 *   is a degeneracy (the dna can combine both ways uniquely).
 * - concA/concB are optional; if both are omitted, the concentration term R·ln(term) is excluded (treated as 0).

 *
 * 			(For hairpins, is this a relavant case?)
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Parameters, Returns, and Errors
 * ──────────────────────────────────────────────────────────────────────────
 * @param {number} sumDeltaH    Total ΔH° in kcal/mol (finite number)
 * @param {number} sumDeltaS    Total ΔS° in cal/K/mol (finite number)
 * @param {number} concA        Initial concentration of strand A, in same units as concB ( > 0 )
 * @param {number} [saltCorrection=0]  Additive entropy correction in cal/K/mol
 * @param {number} [concB]      Initial concentration of strand B ( > 0 if non-self-complementary )
 * @param {boolean} [selfComplementary=false]  Set true for self-complimentary DNA
 *
 * @returns {number} Tm in °C (rounded to TM_DECIMAL_PLACES)
 *
 * @throws {Error}
 *   - If any parameter is missing/invalid
 *   - If any required concentration ≤ 0
 *   - If the concentration term computed for ln(·) is ≤ 0
 *   - If the computed denominator (ΔS° + R·ln(term)) is ~0 or yields non-physical Tm (≤ 0 K)
 */
function calculateTm(
	sumDeltaH,
	sumDeltaS,
	concA,
	concB,
	selfComplementary = false,
	saltCorrection = 0,
) {
	//──────────────────────────────────────────────────────────────────────────//
	//							Parameter Checking								//
	//──────────────────────────────────────────────────────────────────────────//

	// 1) Validate thermodynamic inputs
	if (typeof sumDeltaH !== 'number' || !Number.isFinite(sumDeltaH)) {
		throw new Error(
			`sumDeltaH must be a finite number (kcal/mol). Received: ${sumDeltaH}`,
		);
	}
	if (typeof sumDeltaS !== 'number' || !Number.isFinite(sumDeltaS)) {
		throw new Error(
			`sumDeltaS must be a finite number (cal/K/mol). Received: ${sumDeltaS}`,
		);
	}

	// 2) Concentrations — optional. Validate only if provided.
	const concAProvided = concA !== undefined && concA !== null;
	const concBProvided = concB !== undefined && concB !== null;
	const anyConcProvided = concAProvided || concBProvided;

	if (anyConcProvided) {
		// Self-complementary: only concA is required if any conc is provided
		if (selfComplementary) {
			if (!concAProvided) {
				throw new Error(
					'concA must be provided for self-complementary duplexes when concentrations are supplied.',
				);
			}
			if (
				typeof concA !== 'number' ||
				!Number.isFinite(concA) ||
				concA <= 0
			) {
				throw new Error(
					`concA must be a positive, finite number in µM. Received: ${concA}`,
				);
			}
		} else {
			// Bimolecular: both concA and concB are required if any conc is provided
			if (!concAProvided || !concBProvided) {
				throw new Error(
					'Both concA and concB must be provided for non-self-complementary duplexes when concentrations are supplied.',
				);
			}
			if (
				typeof concA !== 'number' ||
				!Number.isFinite(concA) ||
				concA <= 0
			) {
				throw new Error(
					`concA must be a positive, finite number in µM. Received: ${concA}`,
				);
			}
			if (
				typeof concB !== 'number' ||
				!Number.isFinite(concB) ||
				concB <= 0
			) {
				throw new Error(
					`concB must be a positive, finite number in µM. Received: ${concB}`,
				);
			}
		}
	}

	// 3) Validate boolean flag
	if (typeof selfComplementary !== 'boolean') {
		throw new Error(
			`selfComplementary must be a boolean. Received: ${selfComplementary}`,
		);
	}

	// 4) Validate saltCorrection
	if (
		typeof saltCorrection !== 'number' ||
		!Number.isFinite(saltCorrection)
	) {
		throw new Error(
			`saltCorrection must be a finite number in cal/K/mol. Received: ${saltCorrection}`,
		);
	}

	//──────────────────────────────────────────────────────────────────────────//
	//								Function Logic								//
	//──────────────────────────────────────────────────────────────────────────//

	// Use consistent units: convert ΔH° to cal/mol; ΔS° already in cal/K/mol.
	const DELTA_H_CAL = sumDeltaH * 1000.0;
	const DELTA_S = sumDeltaS;

	// Gas constant in cal/(K·mol)
	const R_CAL = 1.98720425864;

	// Compute the concentration term for the ln() according to the case.
	let lnContribution = 0.0;

	if (anyConcProvided) {
		let concentrationTerm;

		if (selfComplementary) {
			// A + A ↔ AA ⇒ term = Ct = [A] (initial single-strand concentration)
			concentrationTerm = concA * 1e-6; // convert µM → M
		} else {
			// General A + B ↔ AB at 50% duplex of the limiting strand
			const AB50 = Math.min(concA, concB) / 2.0; // µM
			if (AB50 <= 0) {
				throw new Error(
					`Computed [AB]_50 ≤ 0 from inputs (concA=${concA}, concB=${concB}).`,
				);
			}
			const A50 = concA - AB50; // µM
			const B50 = concB - AB50; // µM
			// term = ([A]_50·[B]_50)/[AB]_50  (units µM)
			concentrationTerm = ((A50 * B50) / AB50) * 1e-6; // M
		}

		// The ln() requires a strictly positive argument
		if (!Number.isFinite(concentrationTerm) || concentrationTerm <= 0) {
			throw new Error(
				`Invalid concentration term for ln(): ${concentrationTerm}. ` +
					`Check input concentrations (must yield a positive term).`,
			);
		}

		lnContribution = R_CAL * Math.log(concentrationTerm);
	}

	// Denominator: ΔS° + R·ln(term)
	const denom = DELTA_S + saltCorrection + lnContribution;

	// Avoid singularity / non-physical results
	if (!Number.isFinite(denom)) {
		throw new Error(
			'Non-finite denominator encountered in Tm calculation.',
		);
	}
	if (Math.abs(denom) < 1e-12) {
		throw new Error(
			'Denominator is ~0 (ΔS° + saltCorrection + R·ln(term) ≈ 0), leading to an infinite Tm. ' +
				'Adjust concentrations or thermodynamic parameters.',
		);
	}

	// Tm in Kelvin
	const Tm_K = DELTA_H_CAL / denom;

	if (!Number.isFinite(Tm_K) || Tm_K <= 0) {
		throw new Error(
			`Computed a non-physical Tm (K) = ${Tm_K}. ` +
				`Check signs/magnitudes of ΔH°, ΔS° and concentrations.`,
		);
	}

	// Convert to °C and round
	const Tm_C = Tm_K - 273.15;
	return parseFloat(Tm_C.toFixed(TM_DECIMAL_PLACES));
}


/*****************************************************************************************/
/******************************** Optional Method Helpers *********************************/
/*****************************************************************************************/
/**
 * Rochester hairpin-loop initiation parameters.
 *
 * Computes ΔH° and ΔS° for a hairpin loop of size N using the Rochester
 * parameterization:
 *  - For 3 ≤ N ≤ 30: returns the exact tabulated values from
 *    HAIRPIN_LOOP_PARAMETER_ROCHESTER
 *  - For N > 30: uses the corrected large-N expression:
 *        ΔH°(N) = −14.1  (kcal/mol??)
 *        ΔS°(N) = −[ 14.1 + 6.5 + 0.1*(N − 30) ] * 1000 / 310.15  (cal/mol·K???)
 *
 *
 * Assumptions/Notes
 * -----------------
 *  - Loop size N is an integer count of nucleotides in the loop.
 *  - The 3..30 table contains every integer key
 *
 * @param {number} N
 *   Integer loop size (N ≥ 3). Must be a finite integer.
 *
 * @returns {{ dH: number, dS: number }}
 *   Object with:
 *     - dH: ΔH°(N) in kcal/mol
 *     - dS: ΔS°(N) in cal/(mol·K)
 *
 * @throws {Error}
 *   - If N is not a finite integer.
 *   - If N < 3 (model undefined).
 *
 */
function getRochesterHairpinLoopParams(N) {
	//──────────────────────────────────────────────────────────────────//
	// 						Parameter checking							//
	//──────────────────────────────────────────────────────────────────//
	if (!Number.isFinite(N) || Math.floor(N) !== N) {
		throw new Error('Loop size N must be a finite integer.');
	}
	if (N < 3) {
		throw new Error('Rochester loop initiation is defined for N ≥ 3.');
	}

	//──────────────────────────────────────────────────────────────────//
	// 							Function Logic							//
	//──────────────────────────────────────────────────────────────────//
	// Case 1: Exact table entry for 3..30 (table is dense at unit steps)
	if (N <= 30) {
		const entry = HAIRPIN_LOOP_PARAMETER_ROCHESTER[N];
		// Return a plain object (avoid exposing internal table object)
		return { dH: entry.dH, dS: entry.dS };
	}

	// Case 2: N > 30 → corrected large-N formula
	const dH = -14.1; // kcal/mol
	const dS = -(14.1 + 6.5 + 0.1 * (N - 30)) * (1000 / 310.15); // cal/(mol·K)

	return { dH, dS };
}

/**
 * SantaLucia & Hicks (2004) hairpin-loop initiation parameters.
 * Returns ONLY { dH, dS } with units:
 *   dH in kcal/mol, dS in cal/(mol·K)
 *
 * Rules:
 *  - For exact table sizes (3,4,5,6,7,8,9,10,12,14,16,18,20,25,30):
 *      return stored values (ΔH°=0, ΔS° from table).
 *  - For 3 ≤ N ≤ 30 but not in the table:
 *      linearly interpolate ΔS° between the nearest anchors (ΔH°=0).
 *  - For N > 30:
 *      use corrected asymptotic formula; ΔH°=0,
 *      ΔS°(N) = −[ 6.3 + 1.50*ln(N/30) ] * 1000 / 310.15
 *
 * @param {number} N
 *   Integer loop size (N ≥ 3). Must be a finite integer.
 *
 * @returns {{ dH: number, dS: number }}
 *   Object containing:
 *     - dH: ΔH°(N) in kcal/mol (always 0.0 for Santa Lucia Hicks parameters)
 *     - dS: ΔS°(N) in cal/(mol·K), from table, interpolation, or large-N formula
 *
 * @throws {Error}
 *   - If N is not a finite integer.
 *   - If N < 3 (model is defined only for N ≥ 3).
 */
function getSantaLuciaHicksHairpinParams(N) {
	//──────────────────────────────────────────────────────────────────//
	// 						Parameter checking							//
	//──────────────────────────────────────────────────────────────────//
	if (!Number.isFinite(N) || Math.floor(N) !== N) {
		throw new Error('Loop size N must be a finite integer.');
	}
	if (N < 3) {
		throw new Error(
			'SantaLucia–Hicks loop initiation is defined for N ≥ 3.',
		);
	}

	//──────────────────────────────────────────────────────────────────//
	// 							Function Logic							//
	//──────────────────────────────────────────────────────────────────//

	// Case 1: Exact anchor value present → return directly (ΔH°=0, ΔS° from table)
	if (
		Object.prototype.hasOwnProperty.call(
			HAIRPIN_LOOP_PARAMETERS_SANTA_LUCIA_HICKS,
			N,
		)
	) {
		return HAIRPIN_LOOP_PARAMETERS_SANTA_LUCIA_HICKS[N];
	}

	// Case 2: Interpolate ΔS° for 3 ≤ N ≤ 30 (missing sizes only)
	if (N <= 30) {
		// Parse & sort loop-size keys from the table (integers as strings → numbers)
		const keys = Object.keys(HAIRPIN_LOOP_PARAMETERS_SANTA_LUCIA_HICKS)
			.map(Number)
			.sort((a, b) => a - b);

		// Find nearest lower (lo) and higher (hi) anchors such that lo < N < hi
		let lo = null;
		let hi = null;
		for (let i = 0; i < keys.length - 1; i++) {
			const k0 = keys[i];
			const k1 = keys[i + 1];
			if (k0 < N && N < k1) {
				lo = k0;
				hi = k1;
				break;
			}
		}

		// --- Inline "linear interpolation" on delta S only (delta H=0 across table) ---
		// y(x) = y0 + ((x - x0)/(x1 - x0)) * (y1 - y0)
		const Slo = HAIRPIN_LOOP_PARAMETERS_SANTA_LUCIA_HICKS[lo].dS;
		const Shi = HAIRPIN_LOOP_PARAMETERS_SANTA_LUCIA_HICKS[hi].dS;
		const t = (N - lo) / (hi - lo);
		const dS = Slo + t * (Shi - Slo);

		return { dH: 0.0, dS };
	}

	// Case 3: N > 30
	//   ΔH°(N) = 0
	//   ΔS°(N) = -[ 6.3 + 1.50*ln(N/30) ] * 1000 / 310.15
	// Notes:
	//   - Use natural logarithm (Math.log).
	//   - 310.15 K corresponds to 37 °C.
	const dS = -(6.3 + 1.5 * Math.log(N / 30)) * (1000 / 310.15);

	return { dH: 0.0, dS };
}

/**
 * Returns the Bommarito dangling-end parameters (delta H, delta S) for a given
 * nearest-neighbor and the end on which the dangling base is located (5' or 3').
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Interpretation/Orientation Rules
 * ──────────────────────────────────────────────────────────────────────────
 * - The NN step is written 5'->3' as XY which contains the dangling nucleotide
 * - Orientation rule used:
 *      • 5'-dangling end of form XY -> the dangling base is X
 *      • 3'-dangling end of form XY -> the dangling base is Y
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Assumptions
 * ──────────────────────────────────────────────────────────────────────────
 * - `DANGLING_END_PARAMS` is a deep-frozen object with two
 *   top-level keys: `fivePrime` and `threePrime`, each mapping every possible
 *   dangling end to `{ dH, dS }` in units of kcal/mol and cal/K/mol respectively.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Parameters, Returns, and Errors
 * ──────────────────────────────────────────────────────────────────────────
 * @param {string} step
 *        Two-letter NN step (e.g., "AC"), 5'->3'
 *
 * @param {string} orientation
 *        Accepts '5p'|'3p' or 'fivePrime'|'threePrime' (case-insensitive).
 *
 * @returns {{ dH:number, dS:number }}
 *        Frozen object with ΔH° (kcal/mol) and ΔS° (cal/K/mol).
 *
 * @throws {Error}
 *        If `step` is not exactly two of A/T/C/G, orientation is invalid,
 *        or the pair is not present in the Bommarito table for that orientation.
 */
function getDanglingEndParams(step, orientation) {
	const s = normalizeNNStep(step);
	const o = normalizeDanglingOrientation(orientation);

	const row = DANGLING_END_PARAMS[o][s];

	if (!row) {
		throw new Error(
			`No Bommarito 2000 entry for orientation "${o}" on step "${s}". `,
		);
	}
	return row; // deep-frozen via deepFreeze()
}

/*****************************************************************************************/
/************************************** Small Helpers ************************************/
/*****************************************************************************************/
/**
 * Helper function to freeze all the elements in the table for the loop initiation parameters.
 *
 * @param 		{*} 	obj
 * @returns 	obj		An object with all nested values/objects frozen
 */
function deepFreeze(obj) {
	Object.freeze(obj);
	for (const key of Object.getOwnPropertyNames(obj)) {
		const val = obj[key];
		if (
			val &&
			(typeof val === 'object' || typeof val === 'function') &&
			!Object.isFrozen(val)
		) {
			deepFreeze(val);
		}
	}
	return obj;
}

/**
 * Normalizes an orientation token to the map key: 'fivePrime' | 'threePrime'.
 * @param {string} orientation  '5p'|'3p'|'fivePrime'|'threePrime' (case-insensitive)
 * @returns {'fivePrime'|'threePrime'}
 * @throws {Error} If the token is not recognized.
 */
function normalizeDanglingOrientation(orientation) {
	if (typeof orientation !== 'string') {
		throw new Error(
			`orientation must be a string; got ${typeof orientation}`,
		);
	}
	const o = orientation.trim().toLowerCase();
	if (o === '5p' || o === 'fiveprime') return DANGLING_ORIENTATION.FIVE_PRIME;
	if (o === '3p' || o === 'threeprime')
		return DANGLING_ORIENTATION.THREE_PRIME;
	throw new Error(
		`Unrecognized orientation "${orientation}". Use '5p'|'3p' or 'fivePrime'|'threePrime'.`,
	);
}

/**
 * Ensures the NN query string (for entropy and enthalpy summations) is exactly two DNA bases (A/T/C/G),
 * uppercases it, and returns it.
 *
 * @param {string} step  Two-character nearest-neighbor query string 5'->3'
 * @returns {string}     Uppercased two-letter step
 * @throws {Error}       If not exactly 2 bases A/T/C/G.
 */
function normalizeNNStep(step) {
	if (typeof step !== 'string' || step.length !== 2) {
		throw new Error(
			`NN step must be a 2-character string like "AC"; got: ${step}`,
		);
	}
	const s = step.toUpperCase();
	if (!VALID_BASES.has(s[0]) || !VALID_BASES.has(s[1])) {
		throw new Error(`NN step must contain only A/T/C/G; got: ${step}`);
	}
	return s;
}

/**
 * Builds a canonical "TOP2/BOTTOM2" key after validating each dinucleotide.
 *
 *
 *
 * @param 		{string} 	top2    	Two-letter 5'→3' top-strand dinucleotide
 * @param 		{string} 	bottom2 	Two-letter 3'←5' bottom-strand dinucleotide
 * @returns 	{string} 				e.g., "AC/TG"
 * @throws 		{Error} 				if either part is not exactly two A/T/C/G
 */
function buildTerminalMismatchKey(top2, bottom2) {
	return `${normalizeNNStep(top2)}/${normalizeNNStep(bottom2)}`;
}

/**
 * Parses "TOP2/BOTTOM2" into { top2, bottom2 } using normalizeNNStep.
 * Accepts mixed case and internal whitespace.
 *
 * This function may be helpful in the future.
 *
 * @param 		{string} 			token 		e.g., " aa /  gt "
 * @returns 	{{top2:string, bottom2:string}}
 * @throws 		{Error} 						if malformed or contains invalid bases
 */
function parseTerminalMismatchToken(token) {
	if (typeof token !== 'string') {
		throw new Error(
			`token must be a string like "AA/GT"; got: ${typeof token}`,
		);
	}

	// Remove all types of whitespace and upper case string
	const cleaned = token.replace(/\s+/g, '').toUpperCase();
	const parts = cleaned.split('/');
	if (parts.length !== 2) {
		throw new Error(
			`Malformed terminal mismatch token; expected "NN/NN", got: "${token}"`,
		);
	}
	return {
		top2: normalizeNNStep(parts[0]),
		bottom2: normalizeNNStep(parts[1]),
	};
}

/**
 * Returns the terminal-mismatch delta_H and delta_S for a given NN
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Interpretation
 * ──────────────────────────────────────────────────────────────────────────
 * - Keys represent the two terminal nucleotide paires at the duplex end:
 *     "TOP2/BOTTOM2" where TOP2 is 5'3' on the top strand,
 *     BOTTOM2 is 3'5' on the bottom strand.
 * - Top/bottom is really just a distinction between 5' end mismatches and
 *   3' end mismatches
 * - Both degenerate forms of ("XY/ST", "TS/YX")  are present in the table.
 *
 * ──────────────────────────────────────────────────────────────────────────
 * Parameters, Returns, and Errors
 * ──────────────────────────────────────────────────────────────────────────
 * @param {string} 			top2      					Two-letter top-strand dinucleotide (5'->3')
 * @param {string} 			bottom2   					Two-letter bottom-strand dinucleotide (3'->5')
 *
 * @returns 				{{ dH:number, dS:number }} 	(kcal/mol, cal/K/mol)
 * @throws 					{Error} 					If the NN pair is not found in the table
 */
function getTerminalMismatchParams(top2, bottom2) {
	const key = buildTerminalMismatchKey(top2, bottom2);
	const row = TERMINAL_MISMATCH_PARAMS[key];
	if (!row) {
		throw new Error(`No DNA terminal-mismatch entry for "${key}"`);
	}
	return row; // frozen via deepFreeze()
}

/**
 * Helper function (that we may eventually use) that accepts a single "TOP2/BOTTOM2" token.
 *
 *
 * @param 		{string} 					token 		e.g., "AT/AA"
 * @returns 	{{ dH:number, dS:number }}
 */
function getTerminalMismatchParamsFromToken(token) {
	const { top2, bottom2 } = parseTerminalMismatchToken(token);
	return getTerminalMismatchParams(top2, bottom2);
}


export {
	calculateSnapbackTmRochester,
	calculateSnapbackTmSantaLucia,
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
};
