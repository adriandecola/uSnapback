/*
 * In-house DNA thermodynamic parameters used by the snapback calculators.
 *
 * Units:
 *   dH: kcal/mol
 *   dS: cal/(mol K)
 *
 * Sources are the tables in SW_Methods_Tm_Tool v.2 edit 16:
 *   - SantaLucia & Hicks (2004) matched nearest neighbours and initiation
 *   - Allawi & SantaLucia (1997-1999) single internal mismatches
 *   - SantaLucia & Hicks (2004) hairpin-loop initiation
 */

const row = (dH, dS) => Object.freeze({ dH, dS });

export const DNA_COMPLEMENT = Object.freeze({
	A: 'T',
	T: 'A',
	C: 'G',
	G: 'C',
});

export const VALID_DNA_BASES = Object.freeze(new Set(['A', 'C', 'G', 'T']));

// SantaLucia & Hicks (2004). The key is the 5'->3' dinucleotide on
// the reference strand; its opposing dinucleotide is the Watson-Crick
// complement written 3'->5'. Reverse-equivalent entries are explicit.
export const MATCHED_NN_PARAMS = Object.freeze({
	AA: row(-7.6, -21.3),
	TT: row(-7.6, -21.3),
	AT: row(-7.2, -20.4),
	TA: row(-7.2, -21.3),
	CA: row(-8.5, -22.7),
	TG: row(-8.5, -22.7),
	GT: row(-8.4, -22.4),
	AC: row(-8.4, -22.4),
	CT: row(-7.8, -21.0),
	AG: row(-7.8, -21.0),
	GA: row(-8.2, -22.2),
	TC: row(-8.2, -22.2),
	CG: row(-10.6, -27.2),
	GC: row(-9.8, -24.4),
	GG: row(-8.0, -19.9),
	CC: row(-8.0, -19.9),
});

export const DUPLEX_INITIATION = row(0.2, -5.7);
export const TERMINAL_AT_PENALTY = row(2.2, 6.9);
export const SYMMETRY_CORRECTION = row(0.0, -1.4);

// Each entry supplies both reverse-equivalent tetrads from the source table.
// Key notation is TOP2/BOTTOM2: top is 5'->3', bottom is aligned 3'->5'.
const INTERNAL_MISMATCH_SOURCE_ROWS = Object.freeze([
	// G/T and T/G mismatches — Allawi & SantaLucia (1997)
	['AG/TT', 'TT/GA', 1.0, 0.9],
	['AT/TG', 'GT/TA', -2.5, -8.3],
	['CG/GT', 'TG/GC', -4.1, -11.7],
	['CT/GG', 'GG/TC', -2.8, -8.0],
	['GG/CT', 'TC/GG', 3.3, 10.4],
	['GT/CG', 'GC/TG', -4.4, -12.3],
	['TG/AT', 'TA/GT', -0.1, -1.7],
	['TT/AG', 'GA/TT', -1.3, -5.3],

	// G/A and A/G mismatches — Allawi & SantaLucia (1998)
	['AA/TG', 'GT/AA', -0.6, -2.3],
	['AG/TA', 'AT/GA', -0.7, -2.3],
	['CA/GG', 'GG/AC', -0.7, -2.3],
	['CG/GA', 'AG/GC', -4.0, -13.2],
	['GA/CG', 'GC/AG', -0.6, -1.0],
	['GG/CA', 'AC/GG', 0.5, 3.2],
	['TA/AG', 'GA/AT', 0.7, 0.7],
	['TG/AA', 'AA/GT', 3.0, 7.4],

	// A/C and C/A mismatches at pH 7 — Allawi & SantaLucia (1998)
	['AA/TC', 'CT/AA', 2.3, 4.6],
	['AC/TA', 'AT/CA', 5.3, 14.6],
	['CA/GC', 'CG/AC', 1.9, 3.7],
	['CC/GA', 'AG/CC', 0.6, -0.6],
	['GA/CC', 'CC/AG', 5.2, 14.2],
	['GC/CA', 'AC/CG', -0.7, -3.8],
	['TA/AC', 'CA/AT', 3.4, 8.0],
	['TC/AA', 'AA/CT', 7.6, 20.2],

	// C/T and T/C mismatches — Allawi & SantaLucia (1998)
	['AC/TT', 'TT/CA', 0.7, 0.2],
	['AT/TC', 'CT/TA', -1.2, -6.2],
	['CC/GT', 'TG/CC', -0.8, -4.5],
	['CT/GC', 'CG/TC', -1.5, -6.1],
	['GC/CT', 'TC/CG', 2.3, 5.4],
	['GT/CC', 'CC/TG', 5.2, 13.5],
	['TC/AT', 'TA/CT', 1.2, 0.7],
	['TT/AC', 'CA/TT', 1.0, 0.7],

	// Like-base mismatches — Allawi & SantaLucia (1999)
	['AA/TA', 'AT/AA', 1.2, 1.7],
	['CA/GA', 'AG/AC', -0.9, -4.2],
	['GA/CA', 'AC/AG', -2.9, -9.8],
	['TA/AA', 'AA/AT', 4.7, 12.9],
	['AC/TC', 'CT/CA', 0.0, -4.4],
	['CC/GC', 'CG/CC', -1.5, -7.2],
	['GC/CC', 'CC/CG', 3.6, 8.9],
	['TC/AC', 'CA/CT', 6.1, 16.4],
	['AG/TG', 'GT/GA', -3.1, -9.5],
	['CG/GG', 'GG/GC', -4.9, -15.3],
	['GG/CG', 'GC/GG', -6.0, -15.8],
	['TG/AG', 'GA/GT', 1.6, 3.6],
	['AT/TT', 'TT/TA', -2.7, -10.8],
	['CT/GT', 'TG/TC', -5.0, -15.8],
	['GT/CT', 'TC/TG', -2.2, -8.4],
	['TT/AT', 'TA/TT', 0.2, -1.5],
]);

function buildInternalMismatchParams() {
	const params = {};
	for (const [keyA, keyB, dH, dS] of INTERNAL_MISMATCH_SOURCE_ROWS) {
		const value = row(dH, dS);
		params[keyA] = value;
		params[keyB] = value;
	}
	return Object.freeze(params);
}

export const INTERNAL_MISMATCH_PARAMS = buildInternalMismatchParams();

// Exact SantaLucia & Hicks (2004) table anchors. dG37 is retained so that
// omitted loop sizes can be interpolated in free-energy space before dS is
// derived. Exact listed sizes use the displayed dS values from the write-up.
export const HAIRPIN_LOOP_ANCHORS = Object.freeze({
	3: Object.freeze({ dG37: 3.5, dH: 0.0, dS: -11.3 }),
	4: Object.freeze({ dG37: 3.5, dH: 0.0, dS: -11.3 }),
	5: Object.freeze({ dG37: 3.3, dH: 0.0, dS: -10.6 }),
	6: Object.freeze({ dG37: 4.0, dH: 0.0, dS: -12.9 }),
	7: Object.freeze({ dG37: 4.2, dH: 0.0, dS: -13.5 }),
	8: Object.freeze({ dG37: 4.3, dH: 0.0, dS: -13.9 }),
	9: Object.freeze({ dG37: 4.5, dH: 0.0, dS: -14.5 }),
	10: Object.freeze({ dG37: 4.6, dH: 0.0, dS: -14.8 }),
	12: Object.freeze({ dG37: 5.0, dH: 0.0, dS: -16.1 }),
	14: Object.freeze({ dG37: 5.1, dH: 0.0, dS: -16.4 }),
	16: Object.freeze({ dG37: 5.3, dH: 0.0, dS: -17.1 }),
	18: Object.freeze({ dG37: 5.5, dH: 0.0, dS: -17.7 }),
	20: Object.freeze({ dG37: 5.7, dH: 0.0, dS: -18.4 }),
	25: Object.freeze({ dG37: 6.1, dH: 0.0, dS: -19.7 }),
	30: Object.freeze({ dG37: 6.3, dH: 0.0, dS: -20.3 }),
});

export const GAS_CONSTANT_CAL = 1.98720425864;
export const TEMPERATURE_37_K = 310.15;
