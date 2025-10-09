// tests/script.test.js

import {
	// Primary function
	createSnapback,

	// Secondary functions
	calculateMeltingTempDifferences,
	useForwardPrimer,
	evaluateSnapbackTailMatchingOptions,
	getStemTm,
	getThermoParams,
	createStem,
	buildFinalSnapback,

	// Helper/logic functions
	snvTooCloseToPrimer,
	buildMismatchSequenceForAPI,
	parseTmFromResponse,
	parseThermoParamsFromResponse,
	calculateTm,
	calculateSnapbackTmWittwer,

	// DNA utility functions
	isValidDNASequence,
	isValidSNVObject,
	isValidMismatchObject,
	complementSequence,
	reverseComplement,
	revCompSNV,
	reverseSequence,
	isSelfComplimentary,

	// Loop parameter functions
	getRochesterHairpinLoopParams,
	getSantaLuciaHicksHairpinParams,

	// Dangling end and terminal mismatch helper functions
	getDanglingEndParams,
	normalizeDanglingOrientation,
	normalizeNNStep,
	buildTerminalMismatchKey,
	parseTerminalMismatchToken,
	getTerminalMismatchParamsFromToken,
	getTerminalMismatchParams,

	// Constants
	SNV_BASE_BUFFER,
	NUCLEOTIDE_COMPLEMENT,
	MIN_LOOP_LEN,
	MIN_PRIMER_LEN,
	END_OF_STEM_NUMBER_OF_STRONG_BASE_MISMATCHES_REQUIRED,
	MAX_AMPLICON_LEN,
	HAIRPIN_LOOP_PARAMETER_ROCHESTER,
	HAIRPIN_LOOP_PARAMETERS_SANTA_LUCIA_HICKS,
	DANGLING_END_PARAMS,
	DANGLING_ORIENTATION,
	TERMINAL_MISMATCH_PARAMS,
} from '../dist/script.js';

/****************************************************************/
/*********************** Primary Function ***********************/
/****************************************************************/

describe('createSnapback()', () => {
	test('returns correct snapback sequence, match flags, and Tms for known input case for rs12248560', async () => {
		const targetSeqStrand =
			'ATATTCAGAATAACTAATGTTTGGAAGTTGTTTTGTTTTGCTAAAACAAAGTTTTAGCAAACGATTTTTTTTTTCAAATTTGTGTCTTCTGTTCTCAAAGCATCTCTGATGTAAGAGATAATGCGCCACGATGGGCATCAGAAGACCTCAGCTCAAATCCCAGTTCTGCCAGCTATGAGCTGTGTGGCACCAACAGGTGTC';

		const snvSite = { index: 100, variantBase: 'T' };
		const primerLen = 20;
		const compPrimerLen = 20;
		const targetSnapMeltTemp = 60;

		const result = await createSnapback(
			targetSeqStrand,
			primerLen,
			compPrimerLen,
			snvSite,
			targetSnapMeltTemp
		);

		expect(result.tailOnForwardPrimer).toBe(false); // reverse primer
		expect(result.matchesWild).toBe(true);
		expect(result.snapbackMeltingTms.wildTm).toBeCloseTo(60.45, 2);
		expect(result.snapbackMeltingTms.variantTm).toBeCloseTo(52.508, 2);
		expect(result.snapbackSeq).toBe(
			'AGTGTTCTCAAAGCATCTCTGATGGTGACACCTGTTGGTGCCACAC'
		);
		expect(result.limitingPrimerSeq).toBe('ATATTCAGAATAACTAATGT');
	}, 30000);

	//---------------------------------------------------------------------------------------------------------------------------------//
	const validSeq =
		'GACACCTGTTGGTGCCACACAGCTCATAGCTGGCAGAACTGGGATTTGAGCTGAGGTCTTCTGATGCCCATCGTGGCGCATTATCTCTTACATCAGAGATGCTTTGAGAACAGAAGACACAAATTTGAAAAAAAAAATCGTTTGCTAAAACTTTGTTTTAGCAAAACAAAACAACTTCCAAACATTAGTTATTCTGAATAT';
	const validSNV = { index: 100, variantBase: 'A' };
	const validPrimerLen = 20;
	const validCompPrimerLen = 20;
	const validTm = 60;

	// Invalid targetSeqStrand
	const badTargets = [
		['null', null],
		['non-string', 1234],
		['invalid chars', 'AXTGCT'],
		['lowercase', 'atcg'],
		['array', ['A', 'T', 'C', 'G']],
	];

	for (const [label, target] of badTargets) {
		test(`throws for invalid targetSeqStrand (${label})`, async () => {
			await expect(
				createSnapback(
					target,
					validPrimerLen,
					validCompPrimerLen,
					validSNV,
					validTm
				)
			).rejects.toThrow(/Invalid DNA sequence/i);
		});
	}

	// Invalid Tm
	const badTms = [
		['null', null],
		['negative', -1],
		['zero', 0],
		['non-number', 'hot'],
		['NaN', NaN],
		['Infinity', Infinity],
	];

	for (const [label, tm] of badTms) {
		test(`throws for invalid targetSnapMeltTemp (${label})`, async () => {
			await expect(
				createSnapback(
					validSeq,
					validPrimerLen,
					validCompPrimerLen,
					validSNV,
					tm
				)
			).rejects.toThrow(/targetSnapMeltTemp/i);
		});
	}

	// Invalid primer lengths
	const badPrimerPairs = [
		['primerLen non-int', 12.5, 20],
		['compPrimerLen non-int', 20, 12.2],
		['primerLen negative', -1, 20],
		['compPrimerLen negative', 20, -2],
		['primerLen too short', MIN_PRIMER_LEN - 1, 20],
		['compPrimerLen too short', 20, MIN_PRIMER_LEN - 1],
	];

	for (const [label, primer, compPrimer] of badPrimerPairs) {
		test(`throws for invalid primer lengths (${label})`, async () => {
			await expect(
				createSnapback(validSeq, primer, compPrimer, validSNV, validTm)
			).rejects.toThrow(/primer/i);
		});
	}

	test('throws if total primer lengths exceed sequence length', async () => {
		const longLen = Math.ceil(validSeq.length / 2);
		await expect(
			createSnapback(validSeq, longLen, longLen, validSNV, validTm)
		).rejects.toThrow(/cannot equal or exceed sequence length/i);
	});

	// Invalid SNV object
	const badSNVs = [
		['null', null],
		['missing index', { variantBase: 'A' }],
		['missing variantBase', { index: 5 }],
		['extra key', { index: 5, variantBase: 'A', extra: 1 }],
		['index negative', { index: -1, variantBase: 'A' }],
		['index not integer', { index: 5.5, variantBase: 'A' }],
		['variantBase invalid', { index: 5, variantBase: 'Z' }],
		['variantBase lowercase', { index: 5, variantBase: 'g' }],
		['variantBase too long', { index: 5, variantBase: 'AG' }],
		['variantBase non-string', { index: 5, variantBase: 42 }],
	];

	for (const [label, snv] of badSNVs) {
		test(`throws for invalid snvSite (${label})`, async () => {
			await expect(
				createSnapback(
					validSeq,
					validPrimerLen,
					validCompPrimerLen,
					snv,
					validTm
				)
			).rejects.toThrow(/snvSite|index|variantBase/i);
		});
	}

	test('throws if snvSite.index exceeds sequence length', async () => {
		const snv = { index: validSeq.length + 1, variantBase: 'A' };
		await expect(
			createSnapback(
				validSeq,
				validPrimerLen,
				validCompPrimerLen,
				snv,
				validTm
			)
		).rejects.toThrow(/exceeds sequence length/i);
	});

	test('throws if SNV is too close to primers', async () => {
		const snv = { index: SNV_BASE_BUFFER - 1, variantBase: 'A' };
		await expect(
			createSnapback(
				validSeq,
				validPrimerLen,
				validCompPrimerLen,
				snv,
				validTm
			)
		).rejects.toThrow(/too close to a primer/i);
	});

	// TESTING THE MAX AMPLICON LENGTH TEST
	const makeDNA = (n) => 'ACGT'.repeat(Math.ceil(n / 4)).slice(0, n);
	const primerLen = MIN_PRIMER_LEN;
	const compPrimerLen = MIN_PRIMER_LEN;
	const targetSnapMeltTemp = 60;

	test('throws if amplicon length exceeds MAX_AMPLICON_LEN by 1', async () => {
		const len = MAX_AMPLICON_LEN + 1;
		const targetSeqStrand = makeDNA(len);
		const snvSite = { index: Math.floor(len / 2), variantBase: 'A' };

		await expect(
			createSnapback(
				targetSeqStrand,
				primerLen,
				compPrimerLen,
				snvSite,
				targetSnapMeltTemp
			)
		).rejects.toThrow(/exceeds maximum allowed/i);
	});

	test('throws if amplicon length exceeds MAX_AMPLICON_LEN by a large margin', async () => {
		const len = MAX_AMPLICON_LEN + 250;
		const targetSeqStrand = makeDNA(len);
		const snvSite = { index: Math.floor(len / 2), variantBase: 'C' };

		await expect(
			createSnapback(
				targetSeqStrand,
				primerLen,
				compPrimerLen,
				snvSite,
				targetSnapMeltTemp
			)
		).rejects.toThrow(/exceeds maximum allowed/i);
	});

	test('does not fail the length check when amplicon length equals MAX_AMPLICON_LEN', async () => {
		const len = MAX_AMPLICON_LEN; // boundary: allowed
		const targetSeqStrand = makeDNA(len);
		const snvSite = { index: Math.floor(len / 2), variantBase: 'G' };

		// We only assert it doesn't fail due to the length guard.
		// The function may still reject later for other reasons (that is OK here).
		try {
			await createSnapback(
				targetSeqStrand,
				primerLen,
				compPrimerLen,
				snvSite,
				targetSnapMeltTemp
			);
			// resolved → certainly not the length guard, test passes
		} catch (err) {
			expect(String(err && err.message ? err.message : err)).not.toMatch(
				/exceeds maximum allowed/i
			);
		}
	}, 30000);

	test('does not fail the length check when amplicon length is under MAX_AMPLICON_LEN', async () => {
		const len = MAX_AMPLICON_LEN - 1;
		const targetSeqStrand = makeDNA(len);
		const snvSite = { index: Math.floor(len / 2), variantBase: 'T' };

		try {
			await createSnapback(
				targetSeqStrand,
				primerLen,
				compPrimerLen,
				snvSite,
				targetSnapMeltTemp
			);
		} catch (err) {
			expect(String(err && err.message ? err.message : err)).not.toMatch(
				/exceeds maximum allowed/i
			);
		}
	});
});

/****************************************************************/
/********************* Secondary Functions **********************/
/****************************************************************/

describe('useForwardPrimer()', () => {
	////////////// Logic Test ///////////////
	test('returns correct strand and snapback base for known Tm-diff scenario', async () => {
		const targetSeqStrand =
			'ATATTCAGAATAACTAATGTTTGGAAGTTGTTTTGTTTTGCTAAAACAAAGTTTTAGCAAACGATTTTTTTTTTCAAATTTGTGTCTTCTGTTCTCAAAGCATCTCTGATGTAAGAGATAATGCGCCACGATGGGCATCAGAAGACCTCAGCTCAAATCCCAGTTCTGCCAGCTATGAGCTGTGTGGCACCAACAGGTGTC';
		const snvSite = { index: 100, variantBase: 'T' };

		const result = await useForwardPrimer(targetSeqStrand, snvSite);

		expect(result).toEqual({
			tailOnForwardPrimer: false,
			bestSnapbackTailBaseAtSNV: 'C',
			snapbackTailMatchesWild: true,
		});
	});

	////////////// Parameter Checking ///////////////
	const validTarget = 'ATCGATCGATCGATCGATCG';
	const validSNV = { index: 10, variantBase: 'A' };

	const badTargets = [
		['null', null],
		['non-string', 1234],
		['invalid chars', 'AXTGCT'],
		['non-uppercase', 'atcg'],
		['array instead of string', ['A', 'T', 'C', 'G']],
		['object instead of string', { seq: 'ATCG' }],
	];

	const badSnvs = [
		['null', null],
		['non-object', 'string'],
		['array', [1, 'A']],
		['missing index', { variantBase: 'A' }],
		['missing variantBase', { index: 10 }],
		['extra field', { index: 10, variantBase: 'A', foo: 42 }],
		['non-integer index', { index: 5.5, variantBase: 'A' }],
		['negative index', { index: -1, variantBase: 'A' }],
		['index out of bounds', { index: 100, variantBase: 'A' }],
		['variantBase not a string', { index: 5, variantBase: 7 }],
		['variantBase too long', { index: 5, variantBase: 'AG' }],
		['variantBase invalid', { index: 5, variantBase: 'Z' }],
	];

	for (const [label, badTarget] of badTargets) {
		test(`throws for invalid targetSeqStrand (${label})`, async () => {
			await expect(useForwardPrimer(badTarget, validSNV)).rejects.toThrow(
				/targetSeqStrand/
			);
		});
	}

	for (const [label, badSnv] of badSnvs) {
		test(`throws for invalid snvSite (${label})`, async () => {
			await expect(useForwardPrimer(validTarget, badSnv)).rejects.toThrow(
				/snvSite|index|variantBase/
			);
		});
	}

	////////////// SNV Too Close to Sequence Edges ///////////////
	const SNV_BUFFER = SNV_BASE_BUFFER;
	const makeSeq = (len) => 'A'.repeat(len);

	test('throws if SNV is too close to start (index < SNV_BASE_BUFFER)', async () => {
		const seq = makeSeq(30);
		for (let i = 0; i < SNV_BUFFER; i++) {
			const snv = { index: i, variantBase: 'A' };
			await expect(useForwardPrimer(seq, snv)).rejects.toThrow(
				/is too close to a sequence end/i
			);
		}
	});

	test('throws if SNV is too close to end (index > length - SNV_BASE_BUFFER - 1)', async () => {
		const seq = makeSeq(30);
		const limit = seq.length - SNV_BUFFER - 1;
		for (let i = limit + 1; i < seq.length; i++) {
			const snv = { index: i, variantBase: 'A' };
			await expect(useForwardPrimer(seq, snv)).rejects.toThrow(
				/is too close to a sequence end/i
			);
		}
	});
});

describe('calculateMeltingTempDifferences() — parameter checking', () => {
	// Helper: make a valid uppercase DNA string
	const makeDNA = (n) => 'ACGT'.repeat(Math.ceil(n / 4)).slice(0, n);

	// Valid baseline fixtures — mutate per test
	const validSeq = makeDNA(120); // lengthy enough for various bounds tests
	const validSNV = { index: 60, variantBase: 'A' };
	const validStem = { start: 40, end: 80 }; // includes the SNV at 60
	const validTailFlag = true;

	// Sanity: a wrapper that calls with provided args (so we can re-use)
	const callFn = (
		seq = validSeq,
		snv = validSNV,
		stem = validStem,
		flag = validTailFlag
	) => calculateMeltingTempDifferences(seq, snv, stem, flag);

	// ────────────────────────────────────────────────────────────────────────
	// 1) targetStrandSeqSnapPrimerRefPoint validation
	// ────────────────────────────────────────────────────────────────────────
	const badTargets = [
		['null', null],
		['non-string (number)', 1234],
		['empty string', ''],
		['invalid chars', 'AXTGCT'],
		['lowercase', 'atcg'],
		['array', ['A', 'T', 'C', 'G']],
		['object', { seq: 'ATCG' }],
	];

	for (const [label, badSeq] of badTargets) {
		test(`throws for invalid targetStrandSeqSnapPrimerRefPoint (${label})`, async () => {
			await expect(
				callFn(badSeq, validSNV, validStem, validTailFlag)
			).rejects.toThrow(
				/Invalid targetStrandSeqSnapPrimerRefPoint|DNA string/i
			);
		});
	}

	// ────────────────────────────────────────────────────────────────────────
	// 2) snvSiteSnapPrimerRefPoint validation (shape via helper + bounds here)
	// ────────────────────────────────────────────────────────────────────────
	const badSNVsShape = [
		['null', null],
		['non-object (string)', 'A'],
		['array', [10, 'A']],
		['missing index', { variantBase: 'A' }],
		['missing variantBase', { index: 5 }],
		['index not integer', { index: 5.5, variantBase: 'A' }],
		['variantBase invalid', { index: 5, variantBase: 'Z' }],
		['variantBase lowercase', { index: 5, variantBase: 'g' }],
		['variantBase too long', { index: 5, variantBase: 'AG' }],
		['variantBase non-string', { index: 5, variantBase: 7 }],
	];

	for (const [label, badSnv] of badSNVsShape) {
		test(`throws for invalid snvSiteSnapPrimerRefPoint (${label})`, async () => {
			await expect(
				callFn(validSeq, badSnv, validStem, validTailFlag)
			).rejects.toThrow(/snvSiteSnapPrimerRefPoint|index|variantBase/i);
		});
	}

	test('throws when snvSiteSnapPrimerRefPoint.index is negative (bounds)', async () => {
		await expect(
			callFn(
				validSeq,
				{ index: -1, variantBase: 'A' },
				validStem,
				validTailFlag
			)
		).rejects.toThrow(/out of bounds|index/i);
	});

	test('throws when snvSiteSnapPrimerRefPoint.index equals sequence length (bounds)', async () => {
		const seq = makeDNA(50);
		const snv = { index: seq.length, variantBase: 'C' }; // == SEQ_LEN → OOB
		const stem = { start: 10, end: 20 };
		await expect(callFn(seq, snv, stem, validTailFlag)).rejects.toThrow(
			/out of bounds/i
		);
	});

	test('throws when snvSiteSnapPrimerRefPoint.index greater than sequence length (bounds)', async () => {
		const seq = makeDNA(50);
		const snv = { index: seq.length + 5, variantBase: 'T' };
		const stem = { start: 10, end: 20 };
		await expect(callFn(seq, snv, stem, validTailFlag)).rejects.toThrow(
			/out of bounds/i
		);
	});

	// ────────────────────────────────────────────────────────────────────────
	// 3) bestStemLoc validation (shape, only keys, integer bounds, non-empty)
	// ────────────────────────────────────────────────────────────────────────
	test('throws when bestStemLoc is null', async () => {
		await expect(
			callFn(validSeq, validSNV, null, validTailFlag)
		).rejects.toThrow(/bestStemLoc must be a non-null object/i);
	});

	test('throws when bestStemLoc is non-object', async () => {
		await expect(
			callFn(validSeq, validSNV, 42, validTailFlag)
		).rejects.toThrow(/bestStemLoc must be a non-null object/i);
	});

	test('throws when bestStemLoc is missing start', async () => {
		await expect(
			callFn(validSeq, validSNV, { end: 80 }, validTailFlag)
		).rejects.toThrow(/missing required key "start"/i);
	});

	test('throws when bestStemLoc is missing end', async () => {
		await expect(
			callFn(validSeq, validSNV, { start: 40 }, validTailFlag)
		).rejects.toThrow(/missing required key "end"/i);
	});

	test('throws when bestStemLoc has extra key', async () => {
		await expect(
			callFn(
				validSeq,
				validSNV,
				{ start: 40, end: 80, foo: 1 },
				validTailFlag
			)
		).rejects.toThrow(/unexpected key "foo"/i);
	});

	test('throws when bestStemLoc.start is not an integer', async () => {
		await expect(
			callFn(validSeq, validSNV, { start: 40.5, end: 80 }, validTailFlag)
		).rejects.toThrow(/bestStemLoc\.start must be an integer/i);
	});

	test('throws when bestStemLoc.end is not an integer', async () => {
		await expect(
			callFn(validSeq, validSNV, { start: 40, end: 80.1 }, validTailFlag)
		).rejects.toThrow(/bestStemLoc\.end must be an integer/i);
	});

	test('throws when bestStemLoc.start is NaN', async () => {
		await expect(
			callFn(validSeq, validSNV, { start: NaN, end: 80 }, validTailFlag)
		).rejects.toThrow(/bestStemLoc\.start must be an integer/i);
	});

	test('throws when bestStemLoc.end is Infinity', async () => {
		await expect(
			callFn(
				validSeq,
				validSNV,
				{ start: 40, end: Infinity },
				validTailFlag
			)
		).rejects.toThrow(/bestStemLoc\.end must be an integer/i);
	});

	test('throws when bestStemLoc.start < 0', async () => {
		await expect(
			callFn(validSeq, validSNV, { start: -1, end: 10 }, validTailFlag)
		).rejects.toThrow(/must be ≥ 0/i);
	});

	test('throws when bestStemLoc.end < 0', async () => {
		await expect(
			callFn(validSeq, validSNV, { start: 0, end: -1 }, validTailFlag)
		).rejects.toThrow(/must be ≥ 0/i);
	});

	test('throws when bestStemLoc.start > bestStemLoc.end', async () => {
		await expect(
			callFn(validSeq, validSNV, { start: 50, end: 40 }, validTailFlag)
		).rejects.toThrow(/cannot be greater than/i);
	});

	test('throws when bestStemLoc.end is out of bounds (== SEQ_LEN)', async () => {
		const seq = makeDNA(30);
		const snv = { index: 15, variantBase: 'G' };
		const stem = { start: 10, end: 30 }; // end == SEQ_LEN → OOB
		await expect(callFn(seq, snv, stem, validTailFlag)).rejects.toThrow(
			/out of bounds/i
		);
	});

	test('throws when bestStemLoc.end is out of bounds (> SEQ_LEN)', async () => {
		const seq = makeDNA(30);
		const snv = { index: 15, variantBase: 'C' };
		const stem = { start: 10, end: 35 };
		await expect(callFn(seq, snv, stem, validTailFlag)).rejects.toThrow(
			/out of bounds/i
		);
	});

	test('throws when bestStemLoc spans zero nucleotides (empty)', async () => {
		// The only way to hit "span at least one" guard beyond start>end is to
		// force a degenerate case. Here, we still expect a throw (redundant guard).
		await expect(
			callFn(validSeq, validSNV, { start: 10, end: 9 }, validTailFlag)
		).rejects.toThrow(
			/span at least one nucleotide|cannot be greater than/i
		);
	});

	test('throws when SNV is outside bestStemLoc (below start)', async () => {
		const snv = { index: 9, variantBase: 'A' };
		const stem = { start: 10, end: 20 };
		await expect(
			callFn(validSeq, snv, stem, validTailFlag)
		).rejects.toThrow(/must lie within bestStemLoc/i);
	});

	test('throws when SNV is outside bestStemLoc (above end)', async () => {
		const snv = { index: 21, variantBase: 'T' };
		const stem = { start: 10, end: 20 };
		await expect(
			callFn(validSeq, snv, stem, validTailFlag)
		).rejects.toThrow(/must lie within bestStemLoc/i);
	});

	// Edge inclusivity checks for SNV within the stem bounds
	test('allows SNV exactly at bestStemLoc.start (no throw)', async () => {
		const seq = makeDNA(50);
		const stem = { start: 10, end: 20 };
		const snv = { index: 10, variantBase: 'C' }; // on boundary, valid
		// We only assert that the "must lie within" check does not trigger.
		try {
			await callFn(seq, snv, stem, validTailFlag);
		} catch (err) {
			expect(String(err?.message || err)).not.toMatch(
				/must lie within bestStemLoc/i
			);
		}
	});

	test('allows SNV exactly at bestStemLoc.end (no throw)', async () => {
		const seq = makeDNA(50);
		const stem = { start: 10, end: 20 };
		const snv = { index: 20, variantBase: 'G' }; // on boundary, valid
		try {
			await callFn(seq, snv, stem, validTailFlag);
		} catch (err) {
			expect(String(err?.message || err)).not.toMatch(
				/must lie within bestStemLoc/i
			);
		}
	});

	// ────────────────────────────────────────────────────────────────────────
	// 4) tailOnForwardPrimer validation
	// ────────────────────────────────────────────────────────────────────────
	const badTailFlags = [
		['null', null],
		['number 0', 0],
		['number 1', 1],
		['string "true"', 'true'],
		['string "false"', 'false'],
		['array', []],
		['object', {}],
	];

	for (const [label, badFlag] of badTailFlags) {
		test(`throws for invalid tailOnForwardPrimer (${label})`, async () => {
			await expect(
				callFn(validSeq, validSNV, validStem, badFlag)
			).rejects.toThrow(/tailOnForwardPrimer must be a boolean/i);
		});
	}

	// ────────────────────────────────────────────────────────────────────────
	// 5) Cross-field relational edge cases
	// ────────────────────────────────────────────────────────────────────────
	test('throws when bestStemLoc.end is valid but SNV index equals SEQ_LEN-1 and stem excludes it', async () => {
		const seq = makeDNA(40);
		const snv = { index: seq.length - 1, variantBase: 'A' };
		const stem = { start: 10, end: 20 }; // SNV not inside stem
		await expect(callFn(seq, snv, stem, validTailFlag)).rejects.toThrow(
			/must lie within bestStemLoc/i
		);
	});

	test('throws when bestStemLoc fits sequence, but SNV index out of bounds for that sequence', async () => {
		const seq = makeDNA(25);
		const snv = { index: 30, variantBase: 'T' }; // OOB w.r.t. sequence
		const stem = { start: 5, end: 10 };
		await expect(callFn(seq, snv, stem, validTailFlag)).rejects.toThrow(
			/out of bounds/i
		);
	});
});

describe('evaluateSnapbackTailMatchingOptions()', () => {
	////////////// Logic Testing ///////////////
	test('returns correct snapback tail base and temperature difference for rs12248560, with an SNV base buffer of 4, on the target sequence strand', async () => {
		const initStem = 'AAAGCATCT';
		const mismatchPos = 4;
		const variantBase = 'T';

		const result = await evaluateSnapbackTailMatchingOptions(
			initStem,
			mismatchPos,
			variantBase
		);

		expect(result).toEqual({
			bestSnapbackTailBaseAtSNV: 'A',
			bestTmDifference: 31.54,
			snapbackTailMatchesWild: false,
		});
	});
	test('returns correct snapback base and temperature difference for rs12248560, with an SNV base buffer of 4, on the complement strand', async () => {
		const initStem = 'AAAGCATCT';
		const mismatchPos = 4;
		const variantBase = 'T';
		// Getting the complements for all variables
		const revCompInitStem = reverseComplement(initStem);
		const compMismatchPos = initStem.length - mismatchPos - 1;
		const compVariantBase = NUCLEOTIDE_COMPLEMENT[variantBase];

		const result = await evaluateSnapbackTailMatchingOptions(
			revCompInitStem,
			compMismatchPos,
			compVariantBase
		);

		expect(result).toEqual({
			bestSnapbackTailBaseAtSNV: 'C',
			bestTmDifference: 39.99, // New Tm tool implies it should be 40.00 but thats because thousandth place is rounded
			snapbackTailMatchesWild: true,
		});
	});

	////////////// Parameter Checking ///////////////
	const validSeq = 'GAAAAGGAG';
	const mismatchPos = 4;
	const wild = 'A';
	const variant = 'G';

	test('throws if initStem is not a string', async () => {
		await expect(() =>
			evaluateSnapbackTailMatchingOptions(12345, mismatchPos, variant)
		).rejects.toThrow(/initStem must be/i);
	});

	test('throws if initStem contains invalid characters', async () => {
		await expect(() =>
			evaluateSnapbackTailMatchingOptions(
				'GAAXXGAG',
				mismatchPos,
				variant
			)
		).rejects.toThrow(/initStem must be/i);
	});

	test('throws if mismatchPos is not a number', async () => {
		await expect(() =>
			evaluateSnapbackTailMatchingOptions(validSeq, 'notNum', variant)
		).rejects.toThrow(/must be an integer/i);
	});

	test('throws if mismatchPos is negative', async () => {
		await expect(() =>
			evaluateSnapbackTailMatchingOptions(validSeq, -1, variant)
		).rejects.toThrow(/must be an integer/i);
	});

	test('throws if mismatchPos >= seq.length', async () => {
		await expect(() =>
			evaluateSnapbackTailMatchingOptions(validSeq, 99, variant)
		).rejects.toThrow(/must be an integer between/i);
	});

	test('throws if variantBase is not a string', async () => {
		await expect(() =>
			evaluateSnapbackTailMatchingOptions(validSeq, mismatchPos, 55)
		).rejects.toThrow(/variantBase must be/i);
	});

	test('throws if variantBase is longer than 1 char', async () => {
		await expect(() =>
			evaluateSnapbackTailMatchingOptions(validSeq, mismatchPos, 'TT')
		).rejects.toThrow(/variantBase must be/i);
	});

	test('throws if variantBase is not A/T/C/G', async () => {
		await expect(() =>
			evaluateSnapbackTailMatchingOptions(validSeq, mismatchPos, 'Z')
		).rejects.toThrow(/variantBase must be/i);
	});
});

describe('getStemTm()', () => {
	// Valid parameters
	test('returns 47.27 for sequence "gaaaaggagtgca" with no mismatch', async () => {
		const result = await getStemTm('GAAAAGGAGTGCA');
		expect(result).toBeCloseTo(47.27, 2);
	});

	test('returns 47.27 for sequence "gaaaaggagtgca" with null mismatch', async () => {
		const result = await getStemTm('GAAAAGGAGTGCA', null);
		expect(result).toBeCloseTo(47.27, 2);
	});

	test('returns 37.54 with valid mismatch', async () => {
		const mismatch = { position: 4, type: 'G' };
		const result = await getStemTm('GAAAAGGAGTGCA', mismatch);
		expect(result).toBeCloseTo(37.54, 2);
	});

	test('does NOT throw if mismatch is null', async () => {
		await expect(getStemTm('GAAAAGGAGTGCA', null)).resolves.not.toThrow();
	});

	// Invalid parameters
	test('throws if sequence is not a string', async () => {
		await expect(() => getStemTm(12345)).rejects.toThrow(
			/invalid dna sequence/i
		);
	});

	test('throws if sequence has invalid characters', async () => {
		await expect(() => getStemTm('GAXXXC')).rejects.toThrow(
			/invalid dna sequence/i
		);
	});

	test('throws if mismatch is not an object', async () => {
		await expect(() =>
			getStemTm('GAAAAGGAGTGC', 'not-an-object')
		).rejects.toThrow(/invalid mismatch object/i);
	});

	test('throws if mismatch is an array', async () => {
		await expect(() =>
			getStemTm('GAAAAGGAGTGC', [{ position: 3, type: 'G' }])
		).rejects.toThrow(/invalid mismatch object/i);
	});

	test('throws if mismatch is missing "position"', async () => {
		await expect(() =>
			getStemTm('GAAAAGGAGTGC', { type: 'G' })
		).rejects.toThrow(/invalid mismatch object/i);
	});

	test('throws if mismatch is missing "type"', async () => {
		await expect(() =>
			getStemTm('GAAAAGGAGTGC', { position: 3 })
		).rejects.toThrow(/invalid mismatch object/i);
	});

	test('throws if mismatch.position is out of bounds', async () => {
		await expect(() =>
			getStemTm('GAAAAGGAGTGC', { position: 100, type: 'G' })
		).rejects.toThrow(/exceeds sequence length/i);
	});

	test('throws if mismatch.position is negative', async () => {
		await expect(() =>
			getStemTm('GAAAAGGAGTGC', { position: -1, type: 'G' })
		).rejects.toThrow(/invalid mismatch object/i);
	});

	test('throws if mismatch.position is not an integer', async () => {
		await expect(() =>
			getStemTm('GAAAAGGAGTGC', { position: 2.5, type: 'G' })
		).rejects.toThrow(/invalid mismatch object/i);
	});

	test('throws if mismatch.type is invalid base', async () => {
		await expect(() =>
			getStemTm('GAAAAGGAGTGC', { position: 4, type: 'X' })
		).rejects.toThrow(/invalid mismatch object/i);
	});

	test('throws if mismatch.type is too long', async () => {
		await expect(() =>
			getStemTm('GAAAAGGAGTGC', { position: 4, type: 'GA' })
		).rejects.toThrow(/invalid mismatch object/i);
	});

	test('throws if mismatch has extra keys', async () => {
		await expect(() =>
			getStemTm('GAAAAGGAGTGC', {
				position: 3,
				type: 'A',
				unexpected: true,
			})
		).rejects.toThrow(/invalid mismatch object/i);
	});
});

describe('createStem() – logic test', () => {
	test('returns correct stem location, Tm values, and snapback base for the given sequence', async () => {
		/* ------------------------------------------------------------------
		     Inputs
		   ------------------------------------------------------------------ */
		const targetSeqStrand =
			'GACACCTGTTGGTGCCACACAGCTCATAGCTGGCAGAACTGGGATTTGAGCTGAGGTCTTCTGATGCCCATCGTGGCGCATTATCTCTTACATCAGAGATGCTTTGAGAACAGAAGACACAAATTTGAAAAAAAAAATCGTTTGCTAAAACTTTGTTTTAGCAAAACAAAACAACTTCCAAACATTAGTTATTCTGAATAT';

		const snvSite = { index: 100, variantBase: 'A' }; // the “G”→“A” SNV
		const primerLens = { primerLen: 20, compPrimerLen: 20 };
		const snapbackTailBaseAtSNV = 'C'; // matches wild (complement of G)
		const matchesWild = true;
		const targetSnapMeltTemp = 60;

		/* ------------------------------------------------------------------
		     Call
		   ------------------------------------------------------------------ */
		const result = await createStem(
			targetSeqStrand,
			snvSite,
			primerLens,
			snapbackTailBaseAtSNV,
			matchesWild,
			targetSnapMeltTemp
		);

		/* ------------------------------------------------------------------
		     Expectations
		   ------------------------------------------------------------------ */
		expect(result.bestStemLoc).toEqual({ start: 90, end: 111 });

		expect(result.meltingTemps.wildTm).toBeCloseTo(60.45, 2);
		expect(result.meltingTemps.variantTm).toBeCloseTo(52.508, 2);
	});

	test('returns correct stem location, Tm values, and snapback base when the stem has to move right', async () => {
		/* ------------------------------------------------------------------
		     Inputs
		   ------------------------------------------------------------------ */
		const targetSeqStrand =
			'GACACCTGTTGGTGCCACACAGCTCATAGCTGGCAGAACTGGGATTTGAGCTGAGGTCTTCTGATGCCCATCGTGGCGCATTATCTCTTACATCAGAGATGCTTTGAGAACAGAAGACACAAATTTGAAAAAAAAAATCGTTTGCTAAAACTTTGTTTTAGCAAAACAAAACAACTTCCAAACATTAGTTATTCTGAATAT';

		const snvSite = { index: 100, variantBase: 'A' }; // the “G”→“A” SNV
		const primerLens = { primerLen: 95, compPrimerLen: 20 };
		const snapbackTailBaseAtSNV = 'C'; // matches wild (complement of G)
		const matchesWild = true;
		const targetSnapMeltTemp = 60;

		/* ------------------------------------------------------------------
		     Call
		   ------------------------------------------------------------------ */
		const result = await createStem(
			targetSeqStrand,
			snvSite,
			primerLens,
			snapbackTailBaseAtSNV,
			matchesWild,
			targetSnapMeltTemp
		);

		/* ------------------------------------------------------------------
		     Expectations
		   ------------------------------------------------------------------ */
		expect(result.bestStemLoc.start).toEqual(95);
		// 1. Ensure stem end is greater than 106
		expect(result.bestStemLoc.end).toBeGreaterThan(106);
		// 2. Ensure melt temp is within 2°C of 60
		expect(result.meltingTemps.wildTm).toBeGreaterThanOrEqual(58);
		expect(result.meltingTemps.wildTm).toBeLessThanOrEqual(62);
	});

	test('returns correct stem location, Tm values, and snapback base when the stem has to move left', async () => {
		/* ------------------------------------------------------------------
		     Inputs
		   ------------------------------------------------------------------ */
		const targetSeqStrand =
			'GACACCTGTTGGTGCCACACAGCTCATAGCTGGCAGAACTGGGATTTGAGCTGAGGTCTTCTGATGCCCATCGTGGCGCATTATCTCTTACATCAGAGATGCTTTGAGAACAGAAGACACAAATTTGAAAAAAAAAATCGTTTGCTAAAACTTTGTTTTAGCAAAACAAAACAACTTCCAAACATTAGTTATTCTGAATAT';

		const snvSite = { index: 100, variantBase: 'A' }; // the “G”→“A” SNV
		const primerLens = { primerLen: 20, compPrimerLen: 95 };
		const snapbackTailBaseAtSNV = 'C'; // matches wild (complement of G)
		const matchesWild = true;
		const targetSnapMeltTemp = 60;

		/* ------------------------------------------------------------------
		     Call
		   ------------------------------------------------------------------ */
		const result = await createStem(
			targetSeqStrand,
			snvSite,
			primerLens,
			snapbackTailBaseAtSNV,
			matchesWild,
			targetSnapMeltTemp
		);

		/* ------------------------------------------------------------------
		     Expectations
		   ------------------------------------------------------------------ */
		expect(result.bestStemLoc.end).toEqual(105);
		// 1. Ensure stem start is less than 94
		expect(result.bestStemLoc.start).toBeLessThan(94);
		// 2. Ensure melt temp is within 2°C of 60
		expect(result.meltingTemps.wildTm).toBeGreaterThanOrEqual(58);
		expect(result.meltingTemps.wildTm).toBeLessThanOrEqual(62);
	});

	test('throws error when it cant reach the min snapback melting temperature', async () => {
		/* ------------------------------------------------------------------
		     Inputs
		   ------------------------------------------------------------------ */
		const targetSeqStrand =
			'GACACCTGTTGGTGCCACACAGCTCATAGCTGGCAGAACTGGGATTTGAGCTGAGGTCTTCTGATGCCCATCGTGGCGCATTATCTCTTACATCAGAGATGCTTTGAGAACAGAAGACACAAATTTGAAAAAAAAAATCGTTTGCTAAAACTTTGTTTTAGCAAAACAAAACAACTTCCAAACATTAGTTATTCTGAATAT';

		const snvSite = { index: 100, variantBase: 'A' }; // the “G”→“A” SNV
		const primerLens = { primerLen: 96, compPrimerLen: 96 };
		const snapbackTailBaseAtSNV = 'C'; // matches wild (complement of G)
		const matchesWild = true;
		const targetSnapMeltTemp = 60;

		/* ------------------------------------------------------------------
		     Call
		   ------------------------------------------------------------------ */
		await expect(
			createStem(
				targetSeqStrand,
				snvSite,
				primerLens,
				snapbackTailBaseAtSNV,
				matchesWild,
				targetSnapMeltTemp
			)
		).rejects.toThrow(/Could not meet minimum snapback melting temp/i);
	});
});

describe('createStem() parameter validation', () => {
	const validSeq = 'ATCGATCGATCGATCGAACGCGCGCGCTCGATCGATCGATCGACTATATA';
	const validSNV = { index: 25, variantBase: 'C' };
	const validPrimerLens = { primerLen: 12, compPrimerLen: 12 };
	const validSnapBase = 'G';
	const validMatchesWild = true;
	const validMelt = 67;

	const baseArgs = {
		targetStrandSeqSnapPrimerRefPoint: validSeq,
		snvSiteSnapPrimerRefPoint: validSNV,
		primerLensSnapPrimerRefPoint: validPrimerLens,
		snapbackBaseAtSNV: validSnapBase,
		matchesWild: validMatchesWild,
		desiredSnapbackMeltTempWildType: validMelt,
	};

	const tests = [
		// snvSiteSnapPrimerRefPoint
		[
			'snvSite is null',
			{ snvSiteSnapPrimerRefPoint: null },
			'snvSiteSnapPrimerRefPoint is invalid',
		],
		[
			'snvSite is missing index',
			{ snvSiteSnapPrimerRefPoint: { variantBase: 'A' } },
			'snvSiteSnapPrimerRefPoint is invalid',
		],
		[
			'snvSite index wrong type',
			{ snvSiteSnapPrimerRefPoint: { index: 'a', variantBase: 'A' } },
			'snvSiteSnapPrimerRefPoint is invalid',
		],
		[
			'snvSite index out of range',
			{ snvSiteSnapPrimerRefPoint: { index: 999, variantBase: 'A' } },
			'exceeds sequence length',
		],
		[
			'snvSite variantBase invalid',
			{ snvSiteSnapPrimerRefPoint: { index: 10, variantBase: 'Z' } },
			'snvSiteSnapPrimerRefPoint is invalid',
		],

		// primerLensSnapPrimerRefPoint
		[
			'primerLens is not object',
			{ primerLensSnapPrimerRefPoint: 'notObj' },
			'primerLensSnapPrimerRefPoint',
		],
		[
			'primerLen missing',
			{ primerLensSnapPrimerRefPoint: { compPrimerLen: 12 } },
			'primerLensSnapPrimerRefPoint must be an object with integer properties',
		],
		[
			'compPrimerLen missing',
			{ primerLensSnapPrimerRefPoint: { primerLen: 12 } },
			'primerLensSnapPrimerRefPoint must be an object with integer properties',
		],
		[
			'primerLen < MIN_PRIMER_LEN',
			{
				primerLensSnapPrimerRefPoint: {
					primerLen: MIN_PRIMER_LEN - 1,
					compPrimerLen: MIN_PRIMER_LEN + 1,
				},
			},
			'primerLen must be an integer',
		],
		[
			'compPrimerLen < MIN_PRIMER_LEN',
			{
				primerLensSnapPrimerRefPoint: {
					primerLen: MIN_PRIMER_LEN + 1,
					compPrimerLen: MIN_PRIMER_LEN - 1,
				},
			},
			'compPrimerLen must be an integer',
		],
		[
			'primerLen not int',
			{
				primerLensSnapPrimerRefPoint: {
					primerLen: 12.5,
					compPrimerLen: 12,
				},
			},
			'primerLen must be an integer',
		],
		[
			'compPrimerLen not int',
			{
				primerLensSnapPrimerRefPoint: {
					primerLen: 12,
					compPrimerLen: 12.5,
				},
			},
			'compPrimerLen must be an integer',
		],
		[
			'primerLen negative',
			{
				primerLensSnapPrimerRefPoint: {
					primerLen: -1,
					compPrimerLen: 12,
				},
			},
			'primerLen must be an integer',
		],
		[
			'compPrimerLen negative',
			{
				primerLensSnapPrimerRefPoint: {
					primerLen: 12,
					compPrimerLen: -1,
				},
			},
			'compPrimerLen must be an integer',
		],

		// snapbackTailBaseAtSNV
		[
			'snapbackTailBaseAtSNV invalid base',
			{ snapbackBaseAtSNV: 'Z' },
			'snapbackTailBaseAtSNV',
		],

		// matchesWild
		['matchesWild not boolean', { matchesWild: 'true' }, 'matchesWild'],

		// desiredSnapbackMeltTempWildType
		[
			'desiredSnapbackMeltTempWildType < 0',
			{ desiredSnapbackMeltTempWildType: -5 },
			'must be a positive, finite number',
		],
		[
			'desiredSnapbackMeltTempWildType not a number',
			{ desiredSnapbackMeltTempWildType: 'hot' },
			'must be a positive, finite number',
		],

		// SNV on primer
		[
			'SNV on primer',
			{
				snvSiteSnapPrimerRefPoint: { index: 20, variantBase: 'A' },
				primerLensSnapPrimerRefPoint: {
					primerLen: 20,
					compPrimerLen: 20,
				},
			},
			'SNV at index',
		],
	];

	for (const [label, overrideArgs, errorPattern] of tests) {
		test(`throws for ${label}`, async () => {
			const args = { ...baseArgs, ...overrideArgs }; // overrides valid arg
			await expect(
				createStem(
					args.targetStrandSeqSnapPrimerRefPoint,
					args.snvSiteSnapPrimerRefPoint,
					args.primerLensSnapPrimerRefPoint,
					args.snapbackBaseAtSNV,
					args.matchesWild,
					args.desiredSnapbackMeltTempWildType
				)
			).rejects.toThrow(new RegExp(errorPattern));
		});
	}
});

describe('buildFinalSnapback', () => {
	test('returns correct snapback sequence for known input case', () => {
		const seq =
			'GACACCTGTTGGTGCCACACAGCTCATAGCTGGCAGAACTGGGATTTGAGCTGAGGTCTTCTGATGCCCATCGTGGCGCATTATCTCTTACATCAGAGATGCTTTGAGAACAGAAGACACAAATTTGAAAAAAAAAATCGTTTGCTAAAACTTTGTTTTAGCAAAACAAAACAACTTCCAAACATTAGTTATTCTGAATAT';

		const primers = {
			primerLen: 20,
			compPrimerLen: 20,
		};

		const snv = {
			index: 100,
			variantBase: 'A',
		};

		const stem = {
			start: 90,
			end: 110,
		};

		const tailBaseAtSNV = 'C';

		const expectedSnapback =
			'GAGTTCTCAAAGCATCTCTGATGGTGACACCTGTTGGTGCCACAC';

		const result = buildFinalSnapback(
			seq,
			snv,
			primers,
			stem,
			tailBaseAtSNV
		);
		expect(result).toBe(expectedSnapback);
	});

	// ─────────────────────────────────────────────────────────────────────────────
	// Baseline valid inputs
	// ─────────────────────────────────────────────────────────────────────────────
	const validSeq =
		'ATGCGTGCTAGCTAGCTAGCTAGCTAGCGTATCGATCGTACGTAGCTAGCTAGCGTATCGATC'; // 70 nt
	const validPrimerLens = {
		primerLen: MIN_PRIMER_LEN,
		compPrimerLen: MIN_PRIMER_LEN,
	};
	const validStemLoc = { start: 30, end: 38 };
	const validSNV = { index: 34, variantBase: 'A' };
	const validSnapBase = 'T';

	const baselineCall = () =>
		buildFinalSnapback(
			validSeq,
			validSNV,
			validPrimerLens,
			validStemLoc,
			validSnapBase
		);

	test('baseline call succeeds', () => {
		expect(baselineCall).not.toThrow();
	});

	const callWith = (overrides = {}) => {
		const {
			seq = validSeq,
			snv = validSNV,
			primerLens = validPrimerLens,
			stemLoc = validStemLoc,
			snapBase = validSnapBase,
		} = overrides;
		return () =>
			buildFinalSnapback(seq, snv, primerLens, stemLoc, snapBase);
	};

	// ─────────────────────────────────────────────────────────────────────────────
	// 1. Sequence
	// ─────────────────────────────────────────────────────────────────────────────
	[
		['null', null],
		['number', 1234],
		['lower-case', 'atcg'],
		['invalid chars', 'ATZ#'],
		['array', ['A', 'T']],
	].forEach(([label, bad]) => {
		test(`throws for invalid sequence (${label})`, () => {
			expect(callWith({ seq: bad })).toThrow(/sequence/i);
		});
	});

	// ─────────────────────────────────────────────────────────────────────────────
	// 2. SNV
	// ─────────────────────────────────────────────────────────────────────────────
	const badSNVs = [
		['null', null],
		['missing keys', {}],
		['extra key', { index: 10, variantBase: 'A', foo: 1 }],
		['index negative', { index: -1, variantBase: 'A' }],
		['index non-integer', { index: 2.5, variantBase: 'A' }],
		['index out of bounds', { index: 999, variantBase: 'A' }],
		['variantBase invalid', { index: 10, variantBase: 'Z' }],
		['variantBase long', { index: 10, variantBase: 'AG' }],
	];
	badSNVs.forEach(([label, snv]) => {
		test(`throws for invalid SNV (${label})`, () => {
			expect(callWith({ snv })).toThrow(/SNV/i);
		});
	});

	// ─────────────────────────────────────────────────────────────────────────────
	// 3. snapbackBaseAtSNV
	// ─────────────────────────────────────────────────────────────────────────────
	test('throws for invalid tailBaseAtSNV', () => {
		expect(callWith({ snapBase: 'Z' })).toThrow(/tailBaseAtSNV/);
		expect(callWith({ snapBase: 'AG' })).toThrow(/tailBaseAtSNV/);
		expect(callWith({ snapBase: 42 })).toThrow(/tailBaseAtSNV/);
		expect(callWith({ snapBase: null })).toThrow(/tailBaseAtSNV/);
	});

	// ─────────────────────────────────────────────────────────────────────────────
	// 4. PrimerLensRefPoint
	// ─────────────────────────────────────────────────────────────────────────────
	const badPrimers = [
		['null', null],
		['missing primerLen', { compPrimerLen: MIN_PRIMER_LEN }],
		['missing compPrimerLen', { primerLen: MIN_PRIMER_LEN }],
		['extra key', { primerLen: 20, compPrimerLen: 20, extra: 1 }],
		[
			'primerLen < MIN',
			{ primerLen: MIN_PRIMER_LEN - 1, compPrimerLen: 20 },
		],
		[
			'compPrimerLen < MIN',
			{ primerLen: 20, compPrimerLen: MIN_PRIMER_LEN - 1 },
		],
		['non-int primerLen', { primerLen: 5.5, compPrimerLen: 20 }],
		['non-int compPrimerLen', { primerLen: 20, compPrimerLen: 7.2 }],
		['negative primerLen', { primerLen: -2, compPrimerLen: 20 }],
		['negative compPrimerLen', { primerLen: 20, compPrimerLen: -3 }],
		['sum too large', { primerLen: 35, compPrimerLen: 40 }], // total = 75, same as seq length
	];
	badPrimers.forEach(([label, primerLens]) => {
		test(`throws for invalid primerLens (${label})`, () => {
			expect(callWith({ primerLens })).toThrow(/primer/i);
		});
	});

	// ─────────────────────────────────────────────────────────────────────────────
	// 5. StemLoc
	// ─────────────────────────────────────────────────────────────────────────────
	const badStemLocs = [
		['not object', 'bad'],
		['missing start', { end: 10 }],
		['missing end', { start: 5 }],
		['extra key', { start: 5, end: 10, x: 1 }],
		['start > end', { start: 20, end: 10 }],
		['start < 0', { start: -1, end: 10 }],
		['end < 0', { start: 2, end: -4 }],
		['non-int start', { start: 4.5, end: 10 }],
		['non-int end', { start: 4, end: 10.1 }],
		['end oob', { start: 5, end: validSeq.length }],
		['start oob', { start: validSeq.length, end: validSeq.length }],
	];
	badStemLocs.forEach(([label, stemLoc]) => {
		test(`throws for invalid stemLoc (${label})`, () => {
			expect(callWith({ stemLoc })).toThrow(/stem/i);
		});
	});

	// ─────────────────────────────────────────────────────────────────────────────
	// 6. Functional correctness check
	// ─────────────────────────────────────────────────────────────────────────────
	test('final primer ends with sequence primer region', () => {
		const snap = baselineCall();
		const expectedEnd = validSeq.slice(0, validPrimerLens.primerLen);
		expect(snap.endsWith(expectedEnd)).toBe(true);
	});
});

/****************************************************************/
/*********************** Helper Functions ***********************/
/****************************************************************/

describe('snvTooCloseToPrimer()', () => {
	//These tests assume the SNV_BASE_BUFFER variable will never be unreasonably huge

	// Valid parameter behavior
	test('returns false if SNV is just far away enough from "left-side primer"', () => {
		expect(snvTooCloseToPrimer(5 + SNV_BASE_BUFFER, 5, 10, 50)).toBe(false);
		expect(snvTooCloseToPrimer(6 + SNV_BASE_BUFFER, 6, 10, 60)).toBe(false);
	});

	test('returns false when SNV is just far away enough from "right-side primer', () => {
		const seqLen = 50;
		const primerLen = 5;
		const compPrimerLen = 5;
		const snvIndex = seqLen - compPrimerLen - SNV_BASE_BUFFER - 1;
		expect(
			snvTooCloseToPrimer(snvIndex, primerLen, compPrimerLen, seqLen)
		).toBe(false);
	});

	test('returns false when SNV is comfortably away from both primers', () => {
		expect(
			snvTooCloseToPrimer(
				SNV_BASE_BUFFER + 20,
				5,
				5,
				SNV_BASE_BUFFER + 50
			)
		).toBe(false);
		expect(
			snvTooCloseToPrimer(
				SNV_BASE_BUFFER + 25,
				10,
				10,
				SNV_BASE_BUFFER + 70
			)
		).toBe(false);
	});

	test('returns true if SNV is too close to left primer', () => {
		expect(snvTooCloseToPrimer(5 + SNV_BASE_BUFFER - 2, 5, 10, 50)).toBe(
			true
		);
		expect(snvTooCloseToPrimer(4, 4, 5, 30)).toBe(true);
	});

	test('returns true if SNV is too close to right primer', () => {
		expect(snvTooCloseToPrimer(50 - SNV_BASE_BUFFER - 5, 5, 5, 50)).toBe(
			true
		);
	});

	test('returns true if SNV is both too close to left and right primer', () => {
		expect(snvTooCloseToPrimer(9, 9, 9, 20)).toBe(true);
	});

	test('returns true when SNV is exactly one base inside lower threshold', () => {
		const snvIndex = 5 + SNV_BASE_BUFFER - 1;
		expect(snvTooCloseToPrimer(snvIndex, 5, 5, 50)).toBe(true);
	});

	test('returns true when SNV is exactly one base outside upper threshold', () => {
		const snvIndex = 50 - 5 - SNV_BASE_BUFFER;
		expect(snvTooCloseToPrimer(snvIndex, 5, 5, 50)).toBe(true);
	});

	// Invalid input tests
	test('throws if any input is not a number', () => {
		expect(() => snvTooCloseToPrimer('1', 5, 5, 50)).toThrow(
			'snvIndex must be a finite number'
		);
		expect(() => snvTooCloseToPrimer(1, '5', 5, 50)).toThrow(
			'primerLen must be a finite number'
		);
		expect(() => snvTooCloseToPrimer(1, 5, null, 50)).toThrow(
			'compPrimerLen must be a finite number'
		);
		expect(() => snvTooCloseToPrimer(1, 5, 5, undefined)).toThrow(
			'seqLen must be a finite number'
		);
	});

	test('throws if any input is not an integer', () => {
		expect(() => snvTooCloseToPrimer(1.5, 5, 5, 50)).toThrow(
			'snvIndex must be an integer'
		);
		expect(() => snvTooCloseToPrimer(1, 5.2, 5, 50)).toThrow(
			'primerLen must be an integer'
		);
		expect(() => snvTooCloseToPrimer(1, 5, 5.9, 50)).toThrow(
			'compPrimerLen must be an integer'
		);
		expect(() => snvTooCloseToPrimer(1, 5, 5, 50.01)).toThrow(
			'seqLen must be an integer'
		);
	});

	test('throws for negative values', () => {
		expect(() => snvTooCloseToPrimer(-1, 5, 5, 50)).toThrow(
			'snvIndex must be non-negative'
		);
		expect(() => snvTooCloseToPrimer(1, -5, 5, 50)).toThrow(
			'primerLen must be non-negative'
		);
		expect(() => snvTooCloseToPrimer(1, 5, -5, 50)).toThrow(
			'compPrimerLen must be non-negative'
		);
		expect(() => snvTooCloseToPrimer(1, 5, 5, -50)).toThrow(
			'seqLen must be non-negative'
		);
	});

	test('throws if snvIndex is equal to seqLen', () => {
		expect(() => snvTooCloseToPrimer(50, 5, 5, 50)).toThrow(
			/out of bounds/
		);
	});

	test('throws if snvIndex > seqLen', () => {
		expect(() => snvTooCloseToPrimer(51, 5, 5, 50)).toThrow();
	});

	test('throws if total primer lengths exceed sequence length', () => {
		expect(() => snvTooCloseToPrimer(10, 25, 25, 50)).toThrow(
			/exceed seqLen/
		);
		expect(() => snvTooCloseToPrimer(10, 30, 30, 60)).toThrow(
			/exceed seqLen/
		);
	});

	// Edge cases (can change to throw errors if I want to implement differently)
	test('returns false when primers are 0-length and SNV is within bounds', () => {
		expect(snvTooCloseToPrimer(SNV_BASE_BUFFER + 3, 0, 0, 20)).toBe(false);
		expect(snvTooCloseToPrimer(SNV_BASE_BUFFER, 0, 0, 10)).toBe(false);
	});

	test('returns true when primers are 0 and SNV is too close to sequence edge', () => {
		expect(snvTooCloseToPrimer(0, 0, 0, 10)).toBe(true);
		expect(snvTooCloseToPrimer(10 - SNV_BASE_BUFFER + 2, 0, 0, 10)).toBe(
			true
		);
	});
});

describe('buildMismatchSequenceForAPI()', () => {
	const shortSeq = 'GTCT'; // 4 bp
	const mediumSeq = 'CCTGAGGTTTTAACGGTTGTATGTAGGTGTATATTATCAG'; // 40 bp
	const longSeq =
		'TACTTTGGGATTACACCACCGTCTTTTCGTATGGAACAAGTCTGCCGCGAAACGTGCTGAGTGTACTAACTCGATCTTTCGAGCAAGGAGGTTTTCCGTCCCTGTCTTTTTGTAATCAAA'; // 120 bp
	const veryLongSeq =
		'GCTCGATATCACACCCTGGCCGCGTTCTTTAATTTCATATGAGGCGACTCGTTGCACCCCGTTAACAGAAATCTGAGTGTTATCCCGTTAAGCTAGCGCGTGTCTGCCCGGACACATTTCGCTCCCGAGCCACATTCCTCAAAATCCAACTCGGGGCGCGCGTTGTACTCCTAGACTCCAATATGTAAATACGTGCGTCA'; // 200 bp

	// ====== Valid cases ======
	test('correctly changes base for short sequence at position 0', () => {
		const mismatch = { position: 0, type: 'A' }; // insert T at index 0
		const result = buildMismatchSequenceForAPI(shortSeq, mismatch);
		expect(result.length).toBe(shortSeq.length);
		expect(result[0]).toBe('T');
		expect(result.slice(1)).toBe(shortSeq.slice(1));
	});

	test('correctly changes base for medium sequence in the middle', () => {
		const mismatch = { position: 20, type: 'T' }; // insert A at index 20
		const result = buildMismatchSequenceForAPI(mediumSeq, mismatch);
		expect(result.length).toBe(mediumSeq.length);
		expect(result[20]).toBe('A');
		expect(result.slice(0, 20)).toBe(mediumSeq.slice(0, 20));
		expect(result.slice(21)).toBe(mediumSeq.slice(21));
	});

	test('correctly changes base for long sequence near end', () => {
		const mismatch = { position: 119, type: 'C' }; // insert G at final index
		const result = buildMismatchSequenceForAPI(longSeq, mismatch);
		expect(result.length).toBe(longSeq.length);
		expect(result[119]).toBe('G');
		expect(result.slice(0, 119)).toBe(longSeq.slice(0, 119));
	});

	test('correctly changes base for very long sequence in middle', () => {
		const mismatch = { position: 100, type: 'G' }; // insert C at index 100
		const result = buildMismatchSequenceForAPI(veryLongSeq, mismatch);
		expect(result.length).toBe(veryLongSeq.length);
		expect(result[100]).toBe('C');
		expect(result.slice(0, 100)).toBe(veryLongSeq.slice(0, 100));
		expect(result.slice(101)).toBe(veryLongSeq.slice(101));
	});

	// ====== Invalid input cases ======
	test('throws on invalid base in sequence', () => {
		expect(() =>
			buildMismatchSequenceForAPI('ATBX', { position: 2, type: 'C' })
		).toThrow(/Invalid DNA sequence/i);
	});

	test('throws for non-string sequence', () => {
		expect(() =>
			buildMismatchSequenceForAPI(null, { position: 0, type: 'A' })
		).toThrow(/Invalid DNA sequence/i);
		expect(() =>
			buildMismatchSequenceForAPI(undefined, { position: 0, type: 'A' })
		).toThrow();
		expect(() =>
			buildMismatchSequenceForAPI(42, { position: 0, type: 'A' })
		).toThrow();
	});

	test('throws for missing mismatch keys', () => {
		expect(() => buildMismatchSequenceForAPI('ATCG', {})).toThrow(
			'Invalid mismatch object'
		);
		expect(() =>
			buildMismatchSequenceForAPI('ATCG', { position: 1 })
		).toThrow('Invalid mismatch object');
		expect(() =>
			buildMismatchSequenceForAPI('ATCG', { type: 'A' })
		).toThrow('Invalid mismatch object');
	});

	test('throws for out-of-bounds position', () => {
		expect(() =>
			buildMismatchSequenceForAPI('ATCG', { position: 4, type: 'A' })
		).toThrow(/exceeds sequence length/i);
	});

	test('throws for invalid mismatch type values', () => {
		expect(() =>
			buildMismatchSequenceForAPI('ATCG', { position: 1, type: 'X' })
		).toThrow(/invalid mismatch/i);
		expect(() =>
			buildMismatchSequenceForAPI('ATCG', { position: 1, type: '' })
		).toThrow();
		expect(() =>
			buildMismatchSequenceForAPI('ATCG', { position: 1, type: 'AA' })
		).toThrow();
		expect(() =>
			buildMismatchSequenceForAPI('ATCG', { position: 1, type: 42 })
		).toThrow();
	});

	test('throws for lowercase mismatch type', () => {
		expect(() =>
			buildMismatchSequenceForAPI('ATCG', { position: 1, type: 'g' })
		).toThrow();
	});
});

describe('parseTmFromResponse()', () => {
	// Valid cases
	const validMatchHTML = `
		<html>
		<head></head>
		<body>
			<seq> gaaaaggagtgca </seq>
			<tm> 47.27 </tm>
		</body>
		</html>
	`;

	const validMismatchHTML = `
		<html>
		<head></head>
		<body>
			<seq> gaaaaggagtgca </seq>
			<tm> 47.04 </tm>
			<mmtm> 37.54 </mmtm>
			<warnings> Can only use Santalucia/Hicks parameter set with mismatches. Those parameters will be used here. </warnings>
		</body>
		</html>
	`;

	test('returns 47.27 when mismatch is not passed', () => {
		const result = parseTmFromResponse(validMatchHTML);
		expect(result).toBe(47.27);
	});

	test('returns 47.27 when mismatch is false', () => {
		const result = parseTmFromResponse(validMatchHTML, false);
		expect(result).toBe(47.27);
	});

	test('returns 37.54 when mismatch is true and <mmtm> is present', () => {
		const result = parseTmFromResponse(validMismatchHTML, true);
		expect(result).toBe(37.54);
	});

	// Invalid cases
	test('returns null when <tm> tag is missing and mismatch is false', () => {
		const html = `
			<html>
			<body>
				<seq> ATCG </seq>
			</body>
			</html>
		`;
		expect(parseTmFromResponse(html, false)).toBeNull();
	});

	test('returns null when <mmtm> tag is missing and mismatch is true', () => {
		const html = `
			<html>
			<body>
				<seq> ATCG </seq>
				<tm> 45.5 </tm>
			</body>
			</html>
		`;
		expect(parseTmFromResponse(html, true)).toBeNull();
	});

	test('returns null when <tm> is not a valid number', () => {
		const html = `
			<html><body><tm>not-a-number</tm></body></html>
		`;
		expect(parseTmFromResponse(html)).toBeNull();
	});

	test('returns null when <mmtm> is not a valid number', () => {
		const html = `
			<html><body><mmtm>NaN</mmtm></body></html>
		`;
		expect(parseTmFromResponse(html, true)).toBeNull();
	});

	test('returns null for empty string input', () => {
		expect(parseTmFromResponse('', true)).toBeNull();
		expect(parseTmFromResponse('', false)).toBeNull();
	});

	// 1. Invalid rawHtml values (should all throw)
	const badRawHtmls = [
		['null', null],
		['undefined', undefined],
		['number', 123],
		['object', { html: '<tm>50</tm>' }],
		['array', ['<tm>50</tm>']],
		['boolean', true],
		['function', () => '<tm>50</tm>'],
	];

	for (const [label, badHtml] of badRawHtmls) {
		test(`throws for invalid rawHtml (${label})`, () => {
			expect(() => parseTmFromResponse(badHtml)).toThrow(/rawHtml/i);
		});
	}

	// 2. Invalid mismatch values (should all throw)
	const badMismatches = [
		['string', 'true'],
		['number', 1],
		['object', { mismatch: true }],
		['array', [true]],
		['function', () => true],
	];

	for (const [label, badMismatch] of badMismatches) {
		test(`throws for invalid mismatch (${label})`, () => {
			const html = '<tm>50</tm>';
			expect(() => parseTmFromResponse(html, badMismatch)).toThrow(
				/mismatch/i
			);
		});
	}
});

/****************************************************************/
/******************** DNA Utility Functions *********************/
/****************************************************************/

describe('isValidDNASequence()', () => {
	test('returns true for valid uppercase DNA sequences', () => {
		expect(isValidDNASequence('A')).toBe(true);
		expect(isValidDNASequence('ATCG')).toBe(true);
		expect(isValidDNASequence('GGGTTTAAAACCC')).toBe(true);
	});

	test('returns true for valid single uppercase characters', () => {
		expect(isValidDNASequence('A')).toBe(true);
		expect(isValidDNASequence('T')).toBe(true);
		expect(isValidDNASequence('C')).toBe(true);
		expect(isValidDNASequence('G')).toBe(true);
	});

	test('returns false for lowercase sequences', () => {
		expect(isValidDNASequence('a')).toBe(false);
		expect(isValidDNASequence('atcg')).toBe(false);
		expect(isValidDNASequence('gattaca')).toBe(false);
	});

	test('returns false for lowercase single characters', () => {
		expect(isValidDNASequence('a')).toBe(false);
		expect(isValidDNASequence('t')).toBe(false);
		expect(isValidDNASequence('c')).toBe(false);
		expect(isValidDNASequence('g')).toBe(false);
	});

	test('returns false for mixed-case sequences', () => {
		expect(isValidDNASequence('AtCg')).toBe(false);
		expect(isValidDNASequence('GgGgTtTt')).toBe(false);
	});

	test('returns false for sequences with numbers', () => {
		expect(isValidDNASequence('AT1G')).toBe(false);
		expect(isValidDNASequence('1234')).toBe(false);
		expect(isValidDNASequence('GATT4CA')).toBe(false);
	});

	test('returns false for sequences with special characters', () => {
		expect(isValidDNASequence('A*T^G')).toBe(false);
		expect(isValidDNASequence('@TCG')).toBe(false);
		expect(isValidDNASequence('ACGT!')).toBe(false);
		expect(isValidDNASequence('ACG$T')).toBe(false);
	});

	test('returns false for whitespace-containing strings', () => {
		expect(isValidDNASequence('AT CG')).toBe(false);
		expect(isValidDNASequence(' ATCG')).toBe(false);
		expect(isValidDNASequence('ATCG ')).toBe(false);
		expect(isValidDNASequence('AT\nCG')).toBe(false);
	});

	test('returns false for empty string (no invalid bases)', () => {
		expect(isValidDNASequence('')).toBe(false);
	});

	test('returns false for completely invalid strings', () => {
		expect(isValidDNASequence('xyz')).toBe(false);
		expect(isValidDNASequence('lmnop')).toBe(false);
		expect(isValidDNASequence('....')).toBe(false);
	});

	test('returns false for non-string input (if passed accidentally)', () => {
		expect(isValidDNASequence(1234)).toBe(false);
		expect(isValidDNASequence(null)).toBe(false);
		expect(isValidDNASequence(undefined)).toBe(false);
		expect(isValidDNASequence(true)).toBe(false);
		expect(isValidDNASequence(['A', 'T', 'C', 'G'])).toBe(false);
	});
});

describe('isValidSNVObject', () => {
	// Valid cases
	it('returns true for a valid SNV object', () => {
		expect(isValidSNVObject({ index: 5, variantBase: 'A' })).toBe(true);
		expect(isValidSNVObject({ index: 0, variantBase: 'T' })).toBe(true);
		expect(isValidSNVObject({ index: 123456, variantBase: 'C' })).toBe(
			true
		);
	});

	// Invalid structure
	it('returns false for null, undefined, or non-object values', () => {
		expect(isValidSNVObject(null)).toBe(false);
		expect(isValidSNVObject(undefined)).toBe(false);
		expect(isValidSNVObject(42)).toBe(false);
		expect(isValidSNVObject('A')).toBe(false);
		expect(isValidSNVObject(['index', 5, 'variantBase', 'G'])).toBe(false);
	});

	// Invalid key cases
	it('returns false for missing keys', () => {
		expect(isValidSNVObject({ index: 5 })).toBe(false);
		expect(isValidSNVObject({ variantBase: 'G' })).toBe(false);
		expect(isValidSNVObject({})).toBe(false);
	});

	it('returns false for extra keys', () => {
		expect(
			isValidSNVObject({ index: 5, variantBase: 'C', note: 'extra' })
		).toBe(false);
		expect(isValidSNVObject({ index: 5, variantBase: 'C', pos: 6 })).toBe(
			false
		);
	});

	// Invalid index
	it('returns false for non-integer or negative index values', () => {
		expect(isValidSNVObject({ index: -1, variantBase: 'A' })).toBe(false);
		expect(isValidSNVObject({ index: 1.5, variantBase: 'A' })).toBe(false);
		expect(isValidSNVObject({ index: '5', variantBase: 'G' })).toBe(false);
		expect(isValidSNVObject({ index: NaN, variantBase: 'T' })).toBe(false);
		expect(isValidSNVObject({ index: Infinity, variantBase: 'T' })).toBe(
			false
		);
	});

	// Invalid variantBase
	it('returns false for invalid variantBase characters', () => {
		expect(isValidSNVObject({ index: 3, variantBase: 'U' })).toBe(false);
		expect(isValidSNVObject({ index: 3, variantBase: 'X' })).toBe(false);
		expect(isValidSNVObject({ index: 3, variantBase: 'AT' })).toBe(false);
		expect(isValidSNVObject({ index: 3, variantBase: '' })).toBe(false);
		expect(isValidSNVObject({ index: 3, variantBase: 5 })).toBe(false);
		expect(isValidSNVObject({ index: 3, variantBase: null })).toBe(false);
	});

	it('returns false for lowercase variantBase (must be uppercase)', () => {
		expect(isValidSNVObject({ index: 3, variantBase: 'a' })).toBe(false);
		expect(isValidSNVObject({ index: 3, variantBase: 'g' })).toBe(false);
	});

	// Validity depends on exact key match
	it('returns false when keys are correct but order is wrong in object literal (JS doesn’t care, but test anyway)', () => {
		const snv = { variantBase: 'G', index: 7 }; // order swapped
		expect(isValidSNVObject(snv)).toBe(true); // should still pass
	});
});

describe('complementSequence()', () => {
	test('returns correct complement for simple valid sequences', () => {
		expect(complementSequence('A')).toBe('T');
		expect(complementSequence('T')).toBe('A');
		expect(complementSequence('C')).toBe('G');
		expect(complementSequence('G')).toBe('C');
		expect(complementSequence('ATCG')).toBe('TAGC');
		expect(complementSequence('GGGTTTAAAACCC')).toBe('CCCAAATTTTGGG');
	});

	test('throws error on lowercase input', () => {
		expect(() => complementSequence('atcg')).toThrow(
			'Invalid DNA sequence'
		);
		expect(() => complementSequence('a')).toThrow('Invalid DNA sequence');
	});

	test('throws error on mixed case', () => {
		expect(() => complementSequence('AtCg')).toThrow(
			'Invalid DNA sequence'
		);
	});

	test('throws error on invalid characters', () => {
		expect(() => complementSequence('AT1G')).toThrow(
			'Invalid DNA sequence'
		);
		expect(() => complementSequence('A*T^G')).toThrow(
			'Invalid DNA sequence'
		);
		expect(() => complementSequence('ACGT!')).toThrow(
			'Invalid DNA sequence'
		);
		expect(() => complementSequence('ACGT ')).toThrow(
			'Invalid DNA sequence'
		);
	});

	test('throws error on whitespace', () => {
		expect(() => complementSequence('A T G')).toThrow(
			'Invalid DNA sequence'
		);
		expect(() => complementSequence('AT\nCG')).toThrow(
			'Invalid DNA sequence'
		);
	});

	test('throws error if input is empty', () => {
		expect(() => complementSequence('')).toThrow('Invalid DNA sequence');
	});

	test('throws error for non-string input', () => {
		expect(() => complementSequence(1234)).toThrow('Invalid DNA sequence');
		expect(() => complementSequence(null)).toThrow('Invalid DNA sequence');
		expect(() => complementSequence(undefined)).toThrow(
			'Invalid DNA sequence'
		);
		expect(() => complementSequence(['A', 'T'])).toThrow(
			'Invalid DNA sequence'
		);
	});
});

describe('reverseComplement()', () => {
	// Simple valid tests
	test('returns correct reverse complement for valid sequences', () => {
		expect(reverseComplement('A')).toBe('T');
		expect(reverseComplement('T')).toBe('A');
		expect(reverseComplement('C')).toBe('G');
		expect(reverseComplement('G')).toBe('C');
		expect(reverseComplement('ATCG')).toBe('CGAT');
		expect(reverseComplement('GGGTTTAAA')).toBe('TTTAAACCC');
		expect(reverseComplement('ATGCG')).toBe('CGCAT');
		expect(reverseComplement('TACGT')).toBe('ACGTA');
		expect(reverseComplement('CCGGA')).toBe('TCCGG');
		expect(reverseComplement('TTTGG')).toBe('CCAAA');
	});

	// Long sequences
	test('handles long valid sequences', () => {
		const input = 'ATCGGGTTAACCGGATTCGAGCTTAAACCGGGTTTAAACC';
		const expected = 'GGTTTAAACCCGGTTTAAGCTCGAATCCGGTTAACCCGAT';
		expect(reverseComplement(input)).toBe(expected);
	});

	// Edge case: empty input
	test('throws error for empty input', () => {
		expect(() => reverseComplement('')).toThrow('Invalid DNA sequence');
	});

	// Case sensitivity (should throw on non-uppercase)
	test('throws error on lowercase or mixed-case inputs', () => {
		expect(() => reverseComplement('atcg')).toThrow('Invalid DNA sequence');
		expect(() => reverseComplement('a')).toThrow('Invalid DNA sequence');
		expect(() => reverseComplement('AtCg')).toThrow('Invalid DNA sequence');
		expect(() => reverseComplement('GgGg')).toThrow('Invalid DNA sequence');
	});

	// Invalid characters
	test('throws error on invalid characters or formatting', () => {
		expect(() => reverseComplement('ATXG')).toThrow('Invalid DNA sequence');
		expect(() => reverseComplement('A T')).toThrow('Invalid DNA sequence');
		expect(() => reverseComplement('ACG$')).toThrow('Invalid DNA sequence');
		expect(() => reverseComplement('1234')).toThrow('Invalid DNA sequence');
	});

	// Non-string input
	test('throws error for non-string input', () => {
		expect(() => reverseComplement(42)).toThrow('Invalid DNA sequence');
		expect(() => reverseComplement(null)).toThrow('Invalid DNA sequence');
		expect(() => reverseComplement(undefined)).toThrow(
			'Invalid DNA sequence'
		);
		expect(() => reverseComplement(['A', 'T', 'C', 'G'])).toThrow(
			'Invalid DNA sequence'
		);
	});

	// Palindrome
	test('reverseComplement(reverseComplement(seq)) === seq', () => {
		const sequences = [
			'A',
			'ATCG',
			'GGGTTTAAA',
			'CCGGAATTC',
			'ATCGGGTTAACCGGATTCGAGCTTAAACCGGGTTTAAACC',
		];
		for (const seq of sequences) {
			expect(reverseComplement(reverseComplement(seq))).toBe(seq);
		}
	});
});

describe('revCompSNV()', () => {
	// Valid SNV inputs — variety of lengths
	test('returns correct reverse complement for SNVs in various length sequences', () => {
		expect(revCompSNV({ index: 0, variantBase: 'A' }, 5)).toEqual({
			index: 4,
			variantBase: 'T',
		});
		expect(revCompSNV({ index: 1, variantBase: 'T' }, 6)).toEqual({
			index: 4,
			variantBase: 'A',
		});
		expect(revCompSNV({ index: 2, variantBase: 'C' }, 7)).toEqual({
			index: 4,
			variantBase: 'G',
		});
		expect(revCompSNV({ index: 3, variantBase: 'G' }, 8)).toEqual({
			index: 4,
			variantBase: 'C',
		});
		expect(revCompSNV({ index: 4, variantBase: 'A' }, 9)).toEqual({
			index: 4,
			variantBase: 'T',
		});
		expect(revCompSNV({ index: 5, variantBase: 'T' }, 11)).toEqual({
			index: 5,
			variantBase: 'A',
		});
		expect(revCompSNV({ index: 6, variantBase: 'C' }, 12)).toEqual({
			index: 5,
			variantBase: 'G',
		});
		expect(revCompSNV({ index: 7, variantBase: 'G' }, 13)).toEqual({
			index: 5,
			variantBase: 'C',
		});
		expect(revCompSNV({ index: 8, variantBase: 'A' }, 14)).toEqual({
			index: 5,
			variantBase: 'T',
		});
		expect(revCompSNV({ index: 9, variantBase: 'T' }, 15)).toEqual({
			index: 5,
			variantBase: 'A',
		});
	});

	// Valid SNV inputs
	test('returns correct reverse complement for valid SNVs in 10-length sequence', () => {
		expect(revCompSNV({ index: 0, variantBase: 'A' }, 10)).toEqual({
			index: 9,
			variantBase: 'T',
		});
		expect(revCompSNV({ index: 1, variantBase: 'A' }, 10)).toEqual({
			index: 8,
			variantBase: 'T',
		});
		expect(revCompSNV({ index: 5, variantBase: 'A' }, 10)).toEqual({
			index: 4,
			variantBase: 'T',
		});
		expect(revCompSNV({ index: 9, variantBase: 'A' }, 10)).toEqual({
			index: 0,
			variantBase: 'T',
		});
		expect(revCompSNV({ index: 0, variantBase: 'T' }, 10)).toEqual({
			index: 9,
			variantBase: 'A',
		});
		expect(revCompSNV({ index: 5, variantBase: 'C' }, 10)).toEqual({
			index: 4,
			variantBase: 'G',
		});
		expect(revCompSNV({ index: 7, variantBase: 'G' }, 10)).toEqual({
			index: 2,
			variantBase: 'C',
		});
	});

	// Edge case: index in the middle of even and odd length seqs
	test('returns correct for center index', () => {
		expect(revCompSNV({ index: 4, variantBase: 'G' }, 9)).toEqual({
			index: 4,
			variantBase: 'C',
		});
		expect(revCompSNV({ index: 5, variantBase: 'C' }, 10)).toEqual({
			index: 4,
			variantBase: 'G',
		});
	});

	// Very short sequences
	test('works on shortest valid sequences (length 1-4)', () => {
		expect(revCompSNV({ index: 0, variantBase: 'A' }, 1)).toEqual({
			index: 0,
			variantBase: 'T',
		});
		expect(revCompSNV({ index: 1, variantBase: 'T' }, 3)).toEqual({
			index: 1,
			variantBase: 'A',
		});
		expect(revCompSNV({ index: 2, variantBase: 'G' }, 4)).toEqual({
			index: 1,
			variantBase: 'C',
		});
	});

	// Long sequence case
	test('works on long sequence', () => {
		expect(revCompSNV({ index: 49, variantBase: 'A' }, 100)).toEqual({
			index: 50,
			variantBase: 'T',
		});
		expect(revCompSNV({ index: 0, variantBase: 'G' }, 100)).toEqual({
			index: 99,
			variantBase: 'C',
		});
		expect(revCompSNV({ index: 99, variantBase: 'C' }, 100)).toEqual({
			index: 0,
			variantBase: 'G',
		});
	});

	// Invalid SNV objects
	test('throws for non-object snvSite input', () => {
		expect(() => revCompSNV(null, 10)).toThrow('Invalid SNV object');
		expect(() => revCompSNV(undefined, 10)).toThrow('Invalid SNV object');
		expect(() => revCompSNV([], 10)).toThrow('Invalid SNV object');
		expect(() => revCompSNV(42, 10)).toThrow('Invalid SNV object');
	});

	// SNV object has missing or extra keys
	test('throws if snvSite is missing keys or has extra keys', () => {
		expect(() => revCompSNV({ index: 5 }, 10)).toThrow(
			'Invalid SNV object'
		);
		expect(() => revCompSNV({ variantBase: 'A' }, 10)).toThrow(
			'Invalid SNV object'
		);
		expect(() =>
			revCompSNV({ index: 5, variantBase: 'A', extra: true }, 10)
		).toThrow('Invalid SNV object');
	});

	// Invalid index
	test('throws for invalid index values', () => {
		expect(() => revCompSNV({ index: -1, variantBase: 'A' }, 10)).toThrow(
			'Invalid SNV object'
		);
		expect(() => revCompSNV({ index: 5.5, variantBase: 'C' }, 10)).toThrow(
			'Invalid SNV object'
		);
		expect(() => revCompSNV({ index: '5', variantBase: 'G' }, 10)).toThrow(
			'Invalid SNV object'
		);
	});

	// Invalid base
	test('throws for invalid variantBase', () => {
		expect(() => revCompSNV({ index: 5, variantBase: 'X' }, 10)).toThrow(
			'Invalid SNV object'
		);
		expect(() => revCompSNV({ index: 5, variantBase: 'g' }, 10)).toThrow(
			'Invalid SNV object'
		);
		expect(() => revCompSNV({ index: 5, variantBase: '' }, 10)).toThrow(
			'Invalid SNV object'
		);
		expect(() => revCompSNV({ index: 5, variantBase: 'AA' }, 10)).toThrow(
			'Invalid SNV object'
		);
		expect(() => revCompSNV({ index: 5, variantBase: 1 }, 10)).toThrow(
			'Invalid SNV object'
		);
	});

	test('throws for index out of range', () => {
		expect(() => revCompSNV({ index: 2, variantBase: 'A' }, 1)).toThrow(
			/is out of bounds/
		);
		expect(() => revCompSNV({ index: 10, variantBase: 'T' }, 10)).toThrow(
			/is out of bounds/
		);
		expect(() => revCompSNV({ index: 5000, variantBase: 'C' }, 38)).toThrow(
			/is out of bounds/
		);
		expect(() => revCompSNV({ index: 5, variantBase: 'G' }, 5)).toThrow(
			/is out of bounds/
		);
	});
});

describe('reverseSequence()', () => {
	////////////// Valid Input Tests ///////////////
	test('returns correct reverse for simple sequence', () => {
		expect(reverseSequence('ATCG')).toBe('GCTA');
		expect(reverseSequence('GATTACA')).toBe('ACATTAG');
		expect(reverseSequence('AAAA')).toBe('AAAA');
		expect(reverseSequence('CGTACGTA')).toBe('ATGCATGC');
	});

	test('returns same base for single character', () => {
		expect(reverseSequence('A')).toBe('A');
		expect(reverseSequence('T')).toBe('T');
		expect(reverseSequence('C')).toBe('C');
		expect(reverseSequence('G')).toBe('G');
	});

	////////////// Parameter Checking ///////////////
	const invalidInputs = [
		['null', null],
		['undefined', undefined],
		['number', 1234],
		['array of chars', ['A', 'T', 'C', 'G']],
		['object', { seq: 'ATCG' }],
		['invalid characters', 'AXTG'],
		['lowercase string', 'atcg'],
		['mixed case', 'aTCG'],
		['whitespace', 'A T C G'],
		['special characters', 'AT#G'],
		['newline', 'AT\nCG'],
		['empty sequence', ''],
	];

	for (const [label, badInput] of invalidInputs) {
		test(`throws for invalid input (${label})`, () => {
			expect(() => reverseSequence(badInput)).toThrow(
				/Invalid DNA sequence/
			);
		});
	}
});

describe('calculateSnapbackTmWittwer()', () => {
	////////////////////////////////////////////////////////
	//////////////////// Function Logic ////////////////////
	////////////////////////////////////////////////////////
	const loopLen = 10;

	test('returns correct Tm with no mismatch for wild-type match snapback for rs12248560 on target strand (given by uVariants)', async () => {
		const stemSeq = 'AAAGCATCT';
		const result = await calculateSnapbackTmWittwer(stemSeq, loopLen);
		expect(result).toBe(41.78);
	});

	test('returns correct Tm with mismatch for wild-type match snapback for rs12248560 on target strand (given by uVariants)', async () => {
		const stemSeq = 'AAAGTATCT';
		const mismatch = { position: 4, type: 'G' };
		const result = await calculateSnapbackTmWittwer(
			stemSeq,
			loopLen,
			mismatch
		);
		expect(result).toBe(20.54);
	});

	test('returns correct Tm with no mismatch for variant-type match snapback for rs12248560 on target strand (given by uVariants)', async () => {
		const stemSeq = 'AAAGTATCT';
		const result = await calculateSnapbackTmWittwer(stemSeq, loopLen);
		expect(result).toBe(34.71);
	});

	test('returns correct Tm with mismatch for variant-type match snapback for rs12248560 on target strand (given by uVariants)', async () => {
		const stemSeq = 'AAAGCATCT';
		const mismatch = { position: 4, type: 'A' };
		const result = await calculateSnapbackTmWittwer(
			stemSeq,
			loopLen,
			mismatch
		);
		expect(result).toBe(8.31);
	});

	test('returns correct Tm with no mismatch for wild-type match snapback for rs12248560 on complement strand (given by uVariants)', async () => {
		const stemSeq = 'AGATGCTTT'; // reverse complement of AAAGCATCT
		const result = await calculateSnapbackTmWittwer(stemSeq, loopLen);
		expect(result).toBe(41.78);
	});

	test('returns correct Tm with mismatch for wild-type match snapback for rs12248560 on complement strand (given by uVariants)', async () => {
		const stemSeq = 'AGATACTTT'; // reverse complement of AAAGCATCT
		const mismatch = { position: 4, type: 'C' };
		const result = await calculateSnapbackTmWittwer(
			stemSeq,
			loopLen,
			mismatch
		);
		expect(result).toBe(8.31);
	});

	test('returns correct Tm with no mismatch for variant-type match snapback for rs12248560 on complement strand (given by uVariants)', async () => {
		const stemSeq = 'AGATACTTT'; // reverse complement of AAAGCATCT
		const result = await calculateSnapbackTmWittwer(stemSeq, loopLen);
		expect(result).toBe(34.71);
	});

	test('returns correct Tm with mismatch for variant-type match snapback for rs12248560 on complement strand (given by uVariants)', async () => {
		const stemSeq = 'AGATGCTTT'; // reverse complement of AAAGCATCT
		const mismatch = { position: 4, type: 'T' };
		const result = await calculateSnapbackTmWittwer(
			stemSeq,
			loopLen,
			mismatch
		);
		expect(result).toBe(20.54);
	});

	////////////////////////////////////////////////////////////
	//////////////////// Parameter Checking ////////////////////
	////////////////////////////////////////////////////////////
	const validStem = 'ATCGATCGATCG';
	const validLoopLen = MIN_LOOP_LEN;
	const validMismatch = { position: 4, type: 'T' };

	const badStems = [
		['null', null],
		['number', 12345],
		['non-DNA string', 'AXCGT'],
		['array', ['A', 'T', 'C', 'G']],
		['object', { seq: 'ATCG' }],
		['lowercase string', 'atcg'],
	];

	const badLoopLens = [
		['null', null],
		['non-number', 'loop'],
		['negative', -1],
		['zero', 0],
		['NaN', NaN],
		['Infinity', Infinity],
		['less than MIN_LOOP_LEN', MIN_LOOP_LEN - 1],
	];

	//////////////////////
	// Valid Input Test //
	//////////////////////
	test('returns number for valid input with mismatch', async () => {
		const result = await calculateSnapbackTmWittwer(
			validStem,
			validLoopLen,
			validMismatch
		);
		expect(typeof result).toBe('number');
		expect(Number.isFinite(result)).toBe(true);
	});

	test('returns number for valid input with out a mismatch', async () => {
		const result = await calculateSnapbackTmWittwer(
			validStem,
			validLoopLen
		);
		expect(typeof result).toBe('number');
		expect(Number.isFinite(result)).toBe(true);
	});

	test('returns number for valid input with a null mismatch', async () => {
		const result = await calculateSnapbackTmWittwer(
			validStem,
			validLoopLen,
			null
		);
		expect(typeof result).toBe('number');
		expect(Number.isFinite(result)).toBe(true);
	});

	////////////////////////
	// Invalid Stem Tests //
	////////////////////////
	for (const [label, badStem] of badStems) {
		test(`throws for invalid stemSeq (${label})`, async () => {
			await expect(
				calculateSnapbackTmWittwer(badStem, validLoopLen, validMismatch)
			).rejects.toThrow(/invalid dna sequence/i);
		});
	}

	///////////////////////////
	// Invalid LoopLen Tests //
	///////////////////////////
	for (const [label, badLoop] of badLoopLens) {
		test(`throws for invalid loopLen (${label})`, async () => {
			await expect(
				calculateSnapbackTmWittwer(validStem, badLoop, validMismatch)
			).rejects.toThrow(/loopLen/i);
		});
	}

	////////////////////////////
	// Invalid Mismatch Tests //
	////////////////////////////

	// ────────────────────────────────
	// 1. Wrong Types (not an object)
	// ────────────────────────────────
	const badMismatchTypes = [
		['non-object', 'mismatch'],
		['array', [2, 'A']],
	];

	// ──────────────────────────────────────────────
	// 2. Invalid Shape (missing or extra keys)
	// ──────────────────────────────────────────────
	const badMismatchShapes = [
		['missing position', { type: 'T' }],
		['missing type', { position: 3 }],
		['extra key', { position: 3, type: 'T', extra: 1 }],
	];

	// ─────────────────────────────────────────
	// 3. Invalid Position Field (not valid int)
	// ─────────────────────────────────────────
	const badMismatchPositions = [
		['position not number', { position: '3', type: 'T' }],
		['position not integer', { position: 3.2, type: 'T' }],
		['position negative', { position: -1, type: 'T' }],
	];

	// ──────────────────────────────────────
	// 4. Out-of-bounds Position (SNV rules)
	// ──────────────────────────────────────
	const badMismatchPositionBounds = [
		['position out of bounds', { position: 20, type: 'T' }],
		[
			'not SNV_BASE_BUFFER from start',
			{ position: SNV_BASE_BUFFER - 1, type: 'A' },
		],
		[
			'not SNV_BASE_BUFFER from end',
			{ position: validStem.length - SNV_BASE_BUFFER, type: 'T' },
		],
	];

	// ──────────────────────────────────
	// 5. Invalid Type Field (bad base)
	// ──────────────────────────────────
	const badMismatchTypesOfType = [
		['type invalid base', { position: 2, type: 'Z' }],
		['type wrong type', { position: 2, type: 3 }],
		['type too long', { position: 2, type: 'AG' }],
	];

	// ─────────────────────────────────────
	// Combine and run each sub-category
	// ─────────────────────────────────────
	const badMismatchGroups = [
		...badMismatchTypes,
		...badMismatchShapes,
		...badMismatchPositions,
		...badMismatchTypesOfType,
	];

	for (const [label, badMismatch] of badMismatchGroups) {
		test(`throws for invalid mismatch (${label})`, async () => {
			await expect(
				calculateSnapbackTmWittwer(validStem, validLoopLen, badMismatch)
			).rejects.toThrow(/mismatch/i);
		});
	}

	for (const [label, badMismatch] of badMismatchPositionBounds) {
		test(`throws for out-of-bounds mismatch (${label})`, async () => {
			await expect(
				calculateSnapbackTmWittwer(validStem, validLoopLen, badMismatch)
			).rejects.toThrow(/position/i);
		});
	}
});

describe('isValidMismatchObject()', () => {
	// ─────────────────────────────────────────────────────────────────────────
	// Valid cases
	// ─────────────────────────────────────────────────────────────────────────
	test('returns true for a valid mismatch object with A', () => {
		expect(isValidMismatchObject({ position: 0, type: 'A' })).toBe(true);
	});

	test('returns true for a valid mismatch object with G', () => {
		expect(isValidMismatchObject({ position: 12, type: 'G' })).toBe(true);
	});

	test('returns true for large valid position', () => {
		expect(isValidMismatchObject({ position: 999999, type: 'C' })).toBe(
			true
		);
	});

	// ─────────────────────────────────────────────────────────────────────────
	// Invalid: wrong types
	// ─────────────────────────────────────────────────────────────────────────
	test('returns false if mismatch is null', () => {
		expect(isValidMismatchObject(null)).toBe(false);
	});

	test('returns false if mismatch is an array', () => {
		expect(isValidMismatchObject([1, 2])).toBe(false);
	});

	test('returns false if mismatch is a string', () => {
		expect(isValidMismatchObject('not an object')).toBe(false);
	});

	test('returns false if mismatch is a number', () => {
		expect(isValidMismatchObject(42)).toBe(false);
	});

	// ─────────────────────────────────────────────────────────────────────────
	// Invalid: extra or missing keys
	// ─────────────────────────────────────────────────────────────────────────
	test('returns false if object is missing "position"', () => {
		expect(isValidMismatchObject({ type: 'A' })).toBe(false);
	});

	test('returns false if object is missing "type"', () => {
		expect(isValidMismatchObject({ position: 0 })).toBe(false);
	});

	test('returns false if object has extra key', () => {
		expect(
			isValidMismatchObject({ position: 0, type: 'T', foo: 'bar' })
		).toBe(false);
	});

	// ─────────────────────────────────────────────────────────────────────────
	// Invalid: invalid position
	// ─────────────────────────────────────────────────────────────────────────
	test('returns false if position is negative', () => {
		expect(isValidMismatchObject({ position: -1, type: 'T' })).toBe(false);
	});

	test('returns false if position is a float', () => {
		expect(isValidMismatchObject({ position: 2.5, type: 'C' })).toBe(false);
	});

	test('returns false if position is a string', () => {
		expect(isValidMismatchObject({ position: '5', type: 'G' })).toBe(false);
	});

	test('returns false if position is NaN', () => {
		expect(isValidMismatchObject({ position: NaN, type: 'G' })).toBe(false);
	});

	// ─────────────────────────────────────────────────────────────────────────
	// Invalid: invalid type base
	// ─────────────────────────────────────────────────────────────────────────
	test('returns false if type is a lowercase base', () => {
		expect(isValidMismatchObject({ position: 0, type: 'a' })).toBe(false);
	});

	test('returns false if type is not a valid base', () => {
		expect(isValidMismatchObject({ position: 0, type: 'B' })).toBe(false);
	});

	test('returns false if type is an empty string', () => {
		expect(isValidMismatchObject({ position: 0, type: '' })).toBe(false);
	});

	test('returns false if type is longer than 1 character', () => {
		expect(isValidMismatchObject({ position: 0, type: 'AG' })).toBe(false);
	});

	test('returns false if type is not a string', () => {
		expect(isValidMismatchObject({ position: 0, type: 3 })).toBe(false);
	});
});

/****************************************************************/
/************** Rochester Hairpin Loop Params Tests *************/
/****************************************************************/

describe('getRochesterHairpinLoopParams()', () => {
	// ----------------------------- Helpers ----------------------------- //
	const T_REF = 310.15;
	const largeNdS = (N) => {
		// dS = -[14.1 + 6.5 + 0.1*(N-30)] * 1000 / 310.15
		const base = 14.1 + 6.5 + 0.1 * (N - 30);
		return -((base * 1000) / T_REF);
	};
	const LARGE_N_STEP = -((0.1 * 1000) / T_REF); // ≈ -0.3224246332 cal/(mol·K) per +1 in N

	// ----------------------- Parameter checking ------------------------ //
	test('throws when N is not a finite integer', () => {
		for (const n of [
			null,
			undefined,
			NaN,
			Infinity,
			-Infinity,
			'10',
			12.5,
			{},
			[],
			() => 5,
		]) {
			expect(() => getRochesterHairpinLoopParams(n)).toThrow(
				/finite integer/i
			);
		}
	});

	test('throws when N < 3', () => {
		for (const n of [2, 1, 0, -1, -100]) {
			expect(() => getRochesterHairpinLoopParams(n)).toThrow(
				/N ≥ 3|N >= 3/i
			);
		}
	});

	// ----------------------- Exact table lookups ----------------------- //
	test('returns exact tabulated values for 3 ≤ N ≤ 30 (spot checks)', () => {
		for (const n of [3, 5, 8, 10, 15, 20, 25, 30]) {
			const { dH, dS } = getRochesterHairpinLoopParams(n);
			expect(dH).toBeCloseTo(HAIRPIN_LOOP_PARAMETER_ROCHESTER[n].dH, 6);
			expect(dS).toBeCloseTo(HAIRPIN_LOOP_PARAMETER_ROCHESTER[n].dS, 6);
		}
	});

	test('matches ALL Rochester table entries from N=3..30 exactly (within rounding)', () => {
		for (let n = 3; n <= 30; n++) {
			const { dH, dS } = getRochesterHairpinLoopParams(n);
			expect(dH).toBeCloseTo(HAIRPIN_LOOP_PARAMETER_ROCHESTER[n].dH, 10);
			expect(dS).toBeCloseTo(HAIRPIN_LOOP_PARAMETER_ROCHESTER[n].dS, 10);
		}
	});

	// -------------------------- Large-N branch ------------------------- //
	// Exact expected values from the large-N formula:
	// N=31: dS ≈ -66.74189908108981
	// N=32: dS ≈ -67.06432371433178
	// N=40: dS ≈ -69.64372078026761
	// N=100: dS ≈ -88.9891987747864

	test('N=31: uses large-N formula with correct dH and dS', () => {
		const out = getRochesterHairpinLoopParams(31);
		expect(out.dH).toBeCloseTo(-14.1, 12);
		expect(out.dS).toBeCloseTo(-66.74189908108981, 10);
	});

	test('N=32: uses large-N formula with correct dS', () => {
		const out = getRochesterHairpinLoopParams(32);
		expect(out.dH).toBeCloseTo(-14.1, 12);
		expect(out.dS).toBeCloseTo(-67.06432371433178, 10);
	});

	test('N=40: uses large-N formula with correct dS', () => {
		const out = getRochesterHairpinLoopParams(40);
		expect(out.dH).toBeCloseTo(-14.1, 12);
		expect(out.dS).toBeCloseTo(-69.64372078026761, 10);
	});

	test('N=100: uses large-N formula with correct dS', () => {
		const out = getRochesterHairpinLoopParams(100);
		expect(out.dH).toBeCloseTo(-14.1, 12);
		expect(out.dS).toBeCloseTo(-88.9891987747864, 10);
	});

	test('Large-N branch has linear step size in dS of ≈ -0.3224246332 per +1 in N', () => {
		const dS31 = getRochesterHairpinLoopParams(31).dS;
		const dS32 = getRochesterHairpinLoopParams(32).dS;
		const dS50 = getRochesterHairpinLoopParams(50).dS;
		const dS51 = getRochesterHairpinLoopParams(51).dS;

		expect(dS32 - dS31).toBeCloseTo(LARGE_N_STEP, 10);
		expect(dS51 - dS50).toBeCloseTo(LARGE_N_STEP, 10);
	});

	test('Across the boundary: N=30 (table) vs N=31 (formula) — values are close and ordered', () => {
		const s30 = getRochesterHairpinLoopParams(30).dS; // -66.4
		const s31 = getRochesterHairpinLoopParams(31).dS; // -66.74189...
		// 1) s31 is more negative than s30
		expect(s31).toBeLessThan(s30);
		// 2) The jump magnitude is modest (≈ 0.3224 based on formula step)
		expect(Math.abs(s31 - s30)).toBeLessThan(0.75);
	});

	test('Returns finite numbers at representative N values', () => {
		for (const n of [3, 10, 20, 30, 31, 40, 60, 100]) {
			const { dH, dS } = getRochesterHairpinLoopParams(n);
			expect(Number.isFinite(dH)).toBe(true);
			expect(Number.isFinite(dS)).toBe(true);
		}
	});
});

/*****************************************************/
/***** SantaLucia & Hicks Hairpin Loop Params Tests **/
/*****************************************************/

describe('getSantaLuciaHicksHairpinParams()', () => {
	// ----------------------------- Helpers ----------------------------- //
	const T_REF = 310.15;

	// Large-N entropy (N > 30): dS = -[ 6.3 + 1.50*ln(N/30) ] * 1000 / 310.15
	const largeNdS = (N) => -((6.3 + 1.5 * Math.log(N / 30)) * 1000) / T_REF;

	// Linear interpolation of dS within 3..30 (for non-anchor N only)
	const interpDS = (N) => {
		const keys = Object.keys(HAIRPIN_LOOP_PARAMETERS_SANTA_LUCIA_HICKS)
			.map(Number)
			.sort((a, b) => a - b);

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
		if (lo == null || hi == null) {
			throw new Error(`No interpolation anchors found for N=${N}`);
		}
		const Slo = HAIRPIN_LOOP_PARAMETERS_SANTA_LUCIA_HICKS[lo].dS;
		const Shi = HAIRPIN_LOOP_PARAMETERS_SANTA_LUCIA_HICKS[hi].dS;
		const t = (N - lo) / (hi - lo);
		return Slo + t * (Shi - Slo);
	};

	const isAnchor = (N) =>
		Object.prototype.hasOwnProperty.call(
			HAIRPIN_LOOP_PARAMETERS_SANTA_LUCIA_HICKS,
			N
		);

	// ----------------------- Parameter checking ------------------------ //
	test('throws when N is not a finite integer', () => {
		for (const n of [
			null,
			undefined,
			NaN,
			Infinity,
			-Infinity,
			'10',
			12.3,
			{},
			[],
			() => 7,
		]) {
			expect(() => getSantaLuciaHicksHairpinParams(n)).toThrow(
				/finite integer/i
			);
		}
	});

	test('throws when N < 3', () => {
		for (const n of [2, 1, 0, -1, -100]) {
			expect(() => getSantaLuciaHicksHairpinParams(n)).toThrow(
				/N ≥ 3|N >= 3/i
			);
		}
	});

	// ----------------------- Exact table lookups ----------------------- //
	test('returns exact tabulated values for anchor sizes (spot checks)', () => {
		for (const n of [3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30]) {
			const expected = HAIRPIN_LOOP_PARAMETERS_SANTA_LUCIA_HICKS[n];
			const { dH, dS } = getSantaLuciaHicksHairpinParams(n);
			expect(dH).toBeCloseTo(expected.dH, 12); // always 0.0
			expect(dS).toBeCloseTo(expected.dS, 12);
		}
	});

	// ----------------------- Interpolation branch (3..30 non-anchors) -- //
	test('interpolates dS correctly for non-anchor N in 3..30 (multiple checks)', () => {
		// Choose several non-anchor loop sizes
		const candidates = [11, 13, 15, 17, 19, 21, 22, 23, 24, 26, 27, 28, 29];
		for (const N of candidates) {
			expect(isAnchor(N)).toBe(false);
			const expected_dS = interpDS(N);
			const out = getSantaLuciaHicksHairpinParams(N);
			expect(out.dH).toBeCloseTo(0.0, 12); // ΔH°=0 by definition
			expect(out.dS).toBeCloseTo(expected_dS, 12);
		}
	});

	test('interpolation is bracketed by adjacent anchors (monotone within segment)', () => {
		// Pick a few segments and test an interior N value in each one
		const segments = [
			[10, 12, 11],
			[12, 14, 13],
			[14, 16, 15],
			[18, 20, 19],
			[25, 30, 27], // larger span
		];
		for (const [lo, hi, mid] of segments) {
			const Slo = HAIRPIN_LOOP_PARAMETERS_SANTA_LUCIA_HICKS[lo].dS;
			const Shi = HAIRPIN_LOOP_PARAMETERS_SANTA_LUCIA_HICKS[hi].dS;
			const Smid = getSantaLuciaHicksHairpinParams(mid).dS;
			// Smid should lie within [min(Slo,Shi), max(Slo,Shi)]
			const minS = Math.min(Slo, Shi);
			const maxS = Math.max(Slo, Shi);
			expect(Smid).toBeGreaterThanOrEqual(minS - 1e-9);
			expect(Smid).toBeLessThanOrEqual(maxS + 1e-9);
		}
	});

	// ----------------------- Large-N branch (N > 30) ------------------- //
	// We compute the exact expected values using the large-N formula.
	test('N=31: correct large-N dH and dS', () => {
		const out = getSantaLuciaHicksHairpinParams(31);
		expect(out.dH).toBeCloseTo(0.0, 12);
		expect(out.dS).toBeCloseTo(largeNdS(31), 12);
	});

	test('N=40: correct large-N dS', () => {
		const out = getSantaLuciaHicksHairpinParams(40);
		expect(out.dH).toBeCloseTo(0.0, 12);
		expect(out.dS).toBeCloseTo(largeNdS(40), 12);
	});

	test('N=100: correct large-N dS', () => {
		const out = getSantaLuciaHicksHairpinParams(100);
		expect(out.dH).toBeCloseTo(0.0, 12);
		expect(out.dS).toBeCloseTo(largeNdS(100), 12);
	});

	test('dS decreases (more negative) with N in large-N regime (monotone check)', () => {
		const s31 = getSantaLuciaHicksHairpinParams(31).dS;
		const s60 = getSantaLuciaHicksHairpinParams(60).dS;
		const s100 = getSantaLuciaHicksHairpinParams(100).dS;
		expect(s60).toBeLessThan(s31);
		expect(s100).toBeLessThan(s60);
	});

	// ----------------------- Boundary behavior ------------------------- //
	test('continuity/ordering across boundary: N=30 (table) vs N=31 (formula)', () => {
		const s30 = getSantaLuciaHicksHairpinParams(30).dS;
		const s31 = getSantaLuciaHicksHairpinParams(31).dS;
		// The model does not enforce continuity; just assert sane ordering and smallish jump.
		// Since large-N includes a log term, expect |Δ| < ~2 cal/(mol·K).
		expect(Math.abs(s31 - s30)).toBeLessThan(2.0);
	});

	// ----------------------- General sanity ---------------------------- //
	test('always returns finite numbers and dH === 0 in all regimes', () => {
		for (const N of [3, 4, 10, 11, 15, 20, 25, 29, 30, 31, 40, 75, 100]) {
			const { dH, dS } = getSantaLuciaHicksHairpinParams(N);
			expect(Number.isFinite(dH)).toBe(true);
			expect(Number.isFinite(dS)).toBe(true);
			expect(dH).toBeCloseTo(0.0, 12);
		}
	});
});

describe('Dangling-end parameter helpers (Bommarito 2000)', () => {
	//---------------------------------------------------------------------//
	// Sanity checks on the table shape & immutability
	//---------------------------------------------------------------------//
	test('DANGLING_END_PARAMS is deeply frozen and has 16 steps per orientation', () => {
		expect(Object.isFrozen(DANGLING_END_PARAMS)).toBe(true);
		expect(Object.isFrozen(DANGLING_END_PARAMS.fivePrime)).toBe(true);
		expect(Object.isFrozen(DANGLING_END_PARAMS.threePrime)).toBe(true);

		const fiveKeys = Object.keys(DANGLING_END_PARAMS.fivePrime);
		const threeKeys = Object.keys(DANGLING_END_PARAMS.threePrime);

		expect(fiveKeys).toHaveLength(16);
		expect(threeKeys).toHaveLength(16);

		// All keys should be valid 2-letter DNA steps and rows should be frozen
		for (const k of fiveKeys) {
			expect(k).toMatch(/^[ATCG]{2}$/);
			expect(Object.isFrozen(DANGLING_END_PARAMS.fivePrime[k])).toBe(
				true
			);
		}
		for (const k of threeKeys) {
			expect(k).toMatch(/^[ATCG]{2}$/);
			expect(Object.isFrozen(DANGLING_END_PARAMS.threePrime[k])).toBe(
				true
			);
		}

		// Both orientations should cover all 16 ordered steps
		const ALL_STEPS = (() => {
			const bases = ['A', 'C', 'G', 'T'];
			const out = new Set();
			for (const x of bases) for (const y of bases) out.add(x + y);
			return out;
		})();
		expect(new Set(fiveKeys)).toEqual(ALL_STEPS);
		expect(new Set(threeKeys)).toEqual(ALL_STEPS);
	});

	//---------------------------------------------------------------------//
	// normalizeNNStep(step)
	//---------------------------------------------------------------------//
	test('normalizeNNStep uppercases and validates exactly two DNA letters', () => {
		expect(normalizeNNStep('ac')).toBe('AC');
		expect(normalizeNNStep('tG')).toBe('TG');
		expect(normalizeNNStep('AA')).toBe('AA');
	});

	const badSteps = [
		['non-string', 123],
		['null', null],
		['undefined', undefined],
		['empty', ''],
		['one char', 'A'],
		['three chars', 'AAA'],
		['invalid char X', 'AX'],
		['digit', 'A1'],
		['whitespace padded', ' AC'], // length !== 2
		['tabbed', 'A\t'], // length !== 2
		['object', {}],
		['array', ['A', 'C']],
	];

	for (const [label, s] of badSteps) {
		test(`normalizeNNStep throws for invalid step (${label})`, () => {
			expect(() => normalizeNNStep(s)).toThrow(/NN step/i);
		});
	}

	//---------------------------------------------------------------------//
	// normalizeDanglingOrientation(orientation)
	//---------------------------------------------------------------------//
	test('normalizeDanglingOrientation accepts 5p/3p and fivePrime/threePrime (case-insensitive, trims)', () => {
		expect(normalizeDanglingOrientation('5p')).toBe(
			DANGLING_ORIENTATION.FIVE_PRIME
		);
		expect(normalizeDanglingOrientation('3p')).toBe(
			DANGLING_ORIENTATION.THREE_PRIME
		);

		expect(normalizeDanglingOrientation('fivePrime')).toBe(
			DANGLING_ORIENTATION.FIVE_PRIME
		);
		expect(normalizeDanglingOrientation('ThreePrime')).toBe(
			DANGLING_ORIENTATION.THREE_PRIME
		);

		expect(normalizeDanglingOrientation('  5P  ')).toBe(
			DANGLING_ORIENTATION.FIVE_PRIME
		);
		expect(normalizeDanglingOrientation('\tthreePRIME\n')).toBe(
			DANGLING_ORIENTATION.THREE_PRIME
		);
	});

	const badOrientations = [
		['number', 5],
		['null', null],
		['empty string', ''],
		['random word', 'north'],
		['partial token', 'five'], // not accepted
		['typo', 'threePrim'],
	];

	for (const [label, o] of badOrientations) {
		test(`normalizeDanglingOrientation throws for invalid orientation (${label})`, () => {
			expect(() => normalizeDanglingOrientation(o)).toThrow(
				/orientation/i
			);
		});
	}

	//---------------------------------------------------------------------//
	// getDanglingEndParams(step, orientation)
	//---------------------------------------------------------------------//
	test('getDanglingEndParams returns the frozen row for known examples (case-insensitive step/orientation)', () => {
		// 5′ on AC
		const r1 = getDanglingEndParams('AC', '5p');
		expect(r1).toEqual({ dH: -6.3, dS: -17.1 });
		expect(Object.isFrozen(r1)).toBe(true);

		// Upper/lower case normalization
		const r1b = getDanglingEndParams('ac', 'FivePrime');
		expect(r1b).toEqual(r1);

		// Identity: should be the same object stored in the table (no cloning)
		expect(r1).toBe(DANGLING_END_PARAMS.fivePrime.AC);

		// 3′ on TA
		const r2 = getDanglingEndParams('TA', '3p');
		expect(r2).toEqual({ dH: -0.7, dS: -0.8 });
		expect(Object.isFrozen(r2)).toBe(true);

		// Another 3′ entry with positive values (sanity for sign handling)
		const r3 = getDanglingEndParams('TC', 'threePrime');
		expect(r3).toEqual({ dH: 4.4, dS: 7.3 });

		// 5′ on TA (different from 3′ on TA)
		const r4 = getDanglingEndParams('TA', 'fivePrime');
		expect(r4).toEqual({ dH: -6.9, dS: -20.0 });

		// 3′ on AC (different sign than 5′ on AC)
		const r5 = getDanglingEndParams('AC', '3p');
		expect(r5).toEqual({ dH: 4.7, dS: 14.2 });
	});

	test('getDanglingEndParams result is immutable (attempted mutation does not change values)', () => {
		const r = getDanglingEndParams('GG', 'threePrime');
		const dH = r.dH;
		const dS = r.dS;

		// Try to mutate (in strict mode this would throw; regardless, values must remain)
		try {
			r.dH = 999;
		} catch (_) {}
		try {
			r.extra = 'x';
		} catch (_) {}

		expect(r.dH).toBe(dH);
		expect(r.dS).toBe(dS);
		expect(r).not.toHaveProperty('extra');
	});

	test('getDanglingEndParams throws on invalid orientation tokens', () => {
		expect(() => getDanglingEndParams('AC', 'west')).toThrow(
			/Unrecognized orientation/i
		);
		expect(() => getDanglingEndParams('AC', '')).toThrow(
			/Unrecognized orientation/i
		);
		expect(() => getDanglingEndParams('AC', null)).toThrow(
			/orientation must be a string/i
		);
	});

	test('getDanglingEndParams throws on invalid NN steps (delegates to normalizer)', () => {
		expect(() => getDanglingEndParams('A', '5p')).toThrow(/NN step/i);
		expect(() => getDanglingEndParams('AX', '5p')).toThrow(/NN step/i);
		expect(() => getDanglingEndParams('', '3p')).toThrow(/NN step/i);
		expect(() => getDanglingEndParams(null, '3p')).toThrow(/NN step/i);
	});

	test('getDanglingEndParams matches table values for ALL steps in both orientations', () => {
		for (const [orientationKey, table] of Object.entries(
			DANGLING_END_PARAMS
		)) {
			// orientationKey is 'fivePrime' or 'threePrime'
			const token = orientationKey === 'fivePrime' ? '5p' : '3p';
			for (const [step, expected] of Object.entries(table)) {
				const got = getDanglingEndParams(step, token);
				expect(got).toEqual(expected);

				// Case-insensitive step should match as well
				const lower = step.toLowerCase();
				const gotLower = getDanglingEndParams(
					lower,
					token.toUpperCase()
				);
				expect(gotLower).toEqual(expected);
			}
		}
	});
});

describe('Terminal mismatch lookups (DNA NNDB, Bommarito 2000)', () => {
	//──────────────────────────────────────────────────────────────────//
	// Table shape & immutability
	//──────────────────────────────────────────────────────────────────//
	test('TERMINAL_MISMATCH_PARAMS is deeply frozen with 96 keys', () => {
		expect(Object.isFrozen(TERMINAL_MISMATCH_PARAMS)).toBe(true);

		const keys = Object.keys(TERMINAL_MISMATCH_PARAMS);
		// 48 physical motifs × 2 printed forms = 96 entries
		expect(keys.length).toBe(96);

		for (const k of keys) {
			expect(k).toMatch(/^[ATCG]{2}\/[ATCG]{2}$/);
			const row = TERMINAL_MISMATCH_PARAMS[k];
			expect(Object.isFrozen(row)).toBe(true);
			expect(typeof row.dH).toBe('number');
			expect(typeof row.dS).toBe('number');
		}
	});

	// Quick audit that both printed forms exist for some representative lines
	test('table includes both printed forms for a few representative motifs', () => {
		const pairs = [
			['AA/TA', 'AT/AA'], // A/T neighbor
			['TC/AT', 'TA/CT'], // T/A neighbor
			['CG/GT', 'TG/GC'], // C/G neighbor (with corrected partner)
			['GG/CT', 'TC/GG'], // G/C neighbor
		];
		for (const [a, b] of pairs) {
			expect(TERMINAL_MISMATCH_PARAMS[a]).toBeDefined();
			expect(TERMINAL_MISMATCH_PARAMS[b]).toBeDefined();
			expect(TERMINAL_MISMATCH_PARAMS[a]).toEqual(
				TERMINAL_MISMATCH_PARAMS[b]
			);
		}
	});

	//──────────────────────────────────────────────────────────────────//
	// Normalizers & key helpers (use your normalizeNNStep internally)
	//──────────────────────────────────────────────────────────────────//
	test('buildTerminalMismatchKey normalizes case & validates dinucleotides', () => {
		expect(buildTerminalMismatchKey('ac', 'tg')).toBe('AC/TG');
		expect(() => buildTerminalMismatchKey('A', 'TG')).toThrow(/NN step/i);
		expect(() => buildTerminalMismatchKey('AC', 'T')).toThrow(/NN step/i);
		expect(() => buildTerminalMismatchKey('AX', 'TG')).toThrow(
			/contain only A\/T\/C\/G/i
		);
	});

	test('parseTerminalMismatchToken parses "TOP2/BOTTOM2" with whitespace and mixed case', () => {
		const { top2, bottom2 } = parseTerminalMismatchToken('  tA /  aA ');
		expect(top2).toBe('TA');
		expect(bottom2).toBe('AA');

		expect(() => parseTerminalMismatchToken('TAAA')).toThrow(/Malformed/);
		expect(() => parseTerminalMismatchToken('TA/AX')).toThrow(
			/contain only A\/T\/C\/G/i
		);
		expect(() => parseTerminalMismatchToken(42)).toThrow(/string/i);
	});

	//──────────────────────────────────────────────────────────────────//
	// Core lookup: known values across all neighbor classes
	//──────────────────────────────────────────────────────────────────//
	test('A/T neighbor: AA/TA (and AT/AA) → dH=4.0, dS=15.2', () => {
		const r1 = getTerminalMismatchParams('AA', 'TA');
		expect(r1).toEqual({ dH: 4.0, dS: 15.2 });

		const r2 = getTerminalMismatchParamsFromToken('AT/AA');
		expect(r2).toEqual({ dH: 4.0, dS: 15.2 });
	});

	test('T/A neighbor: TG/AT (and TA/GT) → dH=-5.8, dS=-17.1', () => {
		const r1 = getTerminalMismatchParamsFromToken('TG/AT');
		expect(r1).toEqual({ dH: -5.8, dS: -17.1 });

		const r2 = getTerminalMismatchParams('TA', 'GT');
		expect(r2).toEqual({ dH: -5.8, dS: -17.1 });
	});

	test('C/G neighbor: CG/GT (and TG/GC) → dH=-5.0, dS=-12.9', () => {
		const r1 = getTerminalMismatchParams('CG', 'GT');
		expect(r1).toEqual({ dH: -5.0, dS: -12.9 });

		const r2 = getTerminalMismatchParamsFromToken('TG/GC');
		expect(r2).toEqual({ dH: -5.0, dS: -12.9 });
	});

	test('G/C neighbor: GG/CT (and TC/GG) → dH=-0.9, dS=-0.3', () => {
		const r1 = getTerminalMismatchParams('GG', 'CT');
		expect(r1).toEqual({ dH: -0.9, dS: -0.3 });

		const r2 = getTerminalMismatchParamsFromToken('TC/GG');
		expect(r2).toEqual({ dH: -0.9, dS: -0.3 });
	});

	//──────── Additional spot checks to cover positive/negative enthalpies ───────//
	test('spot checks: positive enthalpy/entropy (A/T neighbor)', () => {
		expect(getTerminalMismatchParams('AT', 'TT')).toEqual({
			dH: 4.3,
			dS: 15.2,
		});
		expect(getTerminalMismatchParams('AG', 'TT')).toEqual({
			dH: 4.3,
			dS: 15.5,
		});
	});

	test('spot checks: negative enthalpy/entropy (T/A, C/G, G/C neighbors)', () => {
		expect(getTerminalMismatchParams('TC', 'AT')).toEqual({
			dH: -5.8,
			dS: -17.1,
		});
		expect(getTerminalMismatchParams('CT', 'GT')).toEqual({
			dH: -5.0,
			dS: -13.2,
		});
		expect(getTerminalMismatchParams('GT', 'CG')).toEqual({
			dH: -5.1,
			dS: -13.5,
		});
	});

	//──────────────────────────────────────────────────────────────────//
	// Immutability of returned rows
	//──────────────────────────────────────────────────────────────────//
	test('returned rows are frozen (attempted mutation has no effect)', () => {
		const row = getTerminalMismatchParams('TT', 'AT'); // -5.8, -17.7
		expect(Object.isFrozen(row)).toBe(true);
		const { dH, dS } = row;
		try {
			row.dH = 999;
		} catch (_) {}
		try {
			row.extra = 'x';
		} catch (_) {}
		expect(row.dH).toBe(dH);
		expect(row.dS).toBe(dS);
		expect(row).not.toHaveProperty('extra');
	});

	//──────────────────────────────────────────────────────────────────//
	// Error handling
	//──────────────────────────────────────────────────────────────────//
	test('unknown motifs throw a clear error', () => {
		// Valid shape but not in table
		expect(() => getTerminalMismatchParams('AA', 'AA')).toThrow(
			/No DNA terminal-mismatch entry/
		);

		// Correct the known transcription error ("T/GC" isn’t valid)
		expect(() => getTerminalMismatchParamsFromToken('T/GC')).toThrow(
			/Malformed|NN step/i
		);
	});

	test('malformed inputs throw (delegated to normalizeNNStep and parser)', () => {
		expect(() => getTerminalMismatchParams('A', 'TA')).toThrow(/NN step/i);
		expect(() => getTerminalMismatchParams('AA', 'T')).toThrow(/NN step/i);
		expect(() => getTerminalMismatchParamsFromToken('AATA')).toThrow(
			/Malformed/
		);
	});
});

describe('calculateTm()', () => {
	const DH = -100.0; // kcal/mol  (example: exothermic duplex)
	const DS = -300.0; // cal/K/mol (example: negative entropy change)
	const Aeq = 1; // 1 µM
	const Beq = 1; // 1 µM

	//──────────────────────────────────────────────────────────────────//
	// Parameter validation
	//──────────────────────────────────────────────────────────────────//
	const badDH = [
		['null', null],
		['NaN', NaN],
		['infinite', Infinity],
		['string', 'hot'],
		['object', {}],
	];

	for (const [label, v] of badDH) {
		test(`throws for invalid sumDeltaH (${label})`, () => {
			expect(() => calculateTm(v, DS, Aeq, Beq, false)).toThrow(
				/sumDeltaH/i
			);
		});
	}

	const badDS = [
		['null', null],
		['NaN', NaN],
		['infinite', -Infinity],
		['string', 'cold'],
		['object', {}],
	];

	for (const [label, v] of badDS) {
		test(`throws for invalid sumDeltaS (${label})`, () => {
			expect(() => calculateTm(DH, v, Aeq, Beq, false)).toThrow(
				/sumDeltaS/i
			);
		});
	}

	const badConcA = [
		['null', null],
		['zero', 0],
		['negative', -1e-6],
		['NaN', NaN],
		['string', '1e-6'],
	];

	for (const [label, v] of badConcA) {
		test(`throws for invalid concA (${label})`, () => {
			expect(() => calculateTm(DH, DS, v, Beq, false)).toThrow(/concA/i);
		});
	}

	const badConcB = [
		['null', null],
		['zero', 0],
		['negative', -2e-6],
		['NaN', NaN],
		['string', '1e-6'],
	];

	for (const [label, v] of badConcB) {
		test(`throws for invalid concB (${label}) when non-self-complementary`, () => {
			expect(() => calculateTm(DH, DS, Aeq, v, false)).toThrow(/concB/i);
		});
	}

	test('throws for invalid selfComplementary flag type', () => {
		// @ts-expect-error - intentionally wrong type for testing
		expect(() => calculateTm(DH, DS, Aeq, Beq, 'yes')).toThrow(
			/selfComplementary/i
		);
	});

	//──────────────────────────────────────────────────────────────────//
	// Core numeric scenarios (rounded to TM_DECIMAL_PLACES)
	//──────────────────────────────────────────────────────────────────//

	test('equimolar strands (A = B): concentration term reduces to Ct/4', () => {
		// Using DH = -100 kcal/mol, DS = -300 cal/K/mol, A = B = 1 uM
		// Expected ≈ 30.96 °C
		const tm = calculateTm(DH, DS, 1, 1, false);
		expect(tm).toBeCloseTo(30.96, 2);
	});

	test('self-complementary case uses concentration term = [A] (Ct)', () => {
		// Using DH = -100 kcal/mol, DS = -300 cal/K/mol, [A] = 1 uM
		// Expected ≈ 32.24 °C
		const tm = calculateTm(DH, DS, 1, undefined, true);
		expect(tm).toBeCloseTo(32.24, 2);
	});

	test('[A] > [B] (limiting strand half-duplexed): term = [A] - [B]/2', () => {
		// A = 2 uM, B = 0.5 uM  → expected ≈ 33.28 °C
		const tm = calculateTm(DH, DS, 2, 0.5, false);
		expect(tm).toBeCloseTo(33.28, 2);
	});

	test('B > A (A limiting): term reduces to [B] - [A]/2 (µM) and gives expected Tm', () => {
		// Choose A = 0.5 µM, B = 2.0 µM.
		// At Tm: term = B - A/2 = 2.0 - 0.25 = 1.75 µM  → 1.75e-6 M
		// With DH = -100 kcal/mol, DS = -300 cal/K/mol ⇒ Tm ≈ 33.28 °C
		const tm = calculateTm(DH, DS, 0.5, 2.0, false);
		expect(tm).toBeCloseTo(33.28, 2);
	});

	test('B > A (A limiting), different numbers: still uses [B] - [A]/2', () => {
		// A = 0.2 µM, B = 3.0 µM → term = 3.0 - 0.1 = 2.9 µM → 2.9e-6 M
		// Expected Tm ≈ 34.22 °C (computed from ΔH/ΔS + R ln(2.9e-6))
		const tm = calculateTm(DH, DS, 0.2, 3.0, false);
		expect(tm).toBeCloseTo(34.22, 2);
	});

	test('symmetry: swapping labels of A and B yields identical Tm', () => {
		const t1 = calculateTm(DH, DS, 2, 0.5, false);
		const t2 = calculateTm(DH, DS, 0.5, 2, false);
		expect(t1).toBeCloseTo(t2, 12);
	});

	//──────────────────────────────────────────────────────────────────//
	// Rounding behavior (to TM_DECIMAL_PLACES)
	//──────────────────────────────────────────────────────────────────//
	test('result is rounded to TM_DECIMAL_PLACES', () => {
		const tm = calculateTm(DH, DS, 1, 1, false);
		// Exact value ≈ 30.956964761... → with 2 decimals it should be 30.96
		expect(tm).toBe(30.96);
	});

	//──────────────────────────────────────────────────────────────────//
	// Non-physical / singular cases
	//──────────────────────────────────────────────────────────────────//
	test('throws if denominator is ~0 (infinite Tm)', () => {
		// Choose DS = 0 and equimolar A=B=2 M so term = (Ct/4) = 1 ⇒ ln(1)=0 ⇒ denom = 0
		expect(() => calculateTm(-10, 0, 2000000, 2000000, false)).toThrow(
			/Denominator is ~0|infinite Tm/i
		);
	});

	test('throws if computed Tm(K) ≤ 0 (non-physical)', () => {
		// Force positive denominator with negative enthalpy large magnitude → negative T
		// Example: DS = +1 cal/K/mol, term = 10 (ln>0), DH = +1 kcal/mol but then T positive; we need T ≤ 0
		// Instead, pick DH very small negative and denom large positive
		expect(() => calculateTm(-0.000001, 100, 1.0, 1.0, false)).toThrow(
			/non-physical Tm/i
		);
	});

	// ──────────────────────────────────────────────────────────────────//
	// saltCorrection behavior
	// ──────────────────────────────────────────────────────────────────//
	test('throws for invalid saltCorrection type', () => {
		// @ts-expect-error - intentionally wrong type for testing
		expect(() => calculateTm(DH, DS, Aeq, Beq, false, '10')).toThrow(
			/saltCorrection/i
		);
	});

	test('explicit saltCorrection = 0 matches default behavior', () => {
		const tDefault = calculateTm(DH, DS, Aeq, Beq, false);
		const tExplicit = calculateTm(DH, DS, Aeq, Beq, false, 0);
		expect(tDefault).toBeCloseTo(tExplicit, 12);
	});

	test('positive saltCorrection raises Tm (equimolar A=B=1 µM)', () => {
		const tm = calculateTm(DH, DS, 1, 1, false, 10); // +10 cal/K/mol
		expect(tm).toBeCloseTo(40.5, 2);
	});

	test('negative saltCorrection lowers Tm (equimolar A=B=1 µM)', () => {
		const tm = calculateTm(DH, DS, 1, 1, false, -10); // -10 cal/K/mol
		expect(tm).toBeCloseTo(21.98, 2);
	});

	test('self-complementary also applies saltCorrection', () => {
		// [A]=1 µM, saltCorrection = +10 cal/K/mol
		const tm = calculateTm(DH, DS, 1, undefined, true, 10);
		expect(tm).toBeCloseTo(41.86, 2);
	});

	test('saltCorrection can drive denominator ~0 and throws', () => {
		// For A=B=1 µM: term = Ct/4 = 0.5 µM ⇒ ln(5e-7) ≈ -14.508657739
		// Choose saltCorrection ≈ -DS - R*ln(term) to cancel denom.
		const sCloseToZero = -DS - 1.98720425864 * Math.log(0.5e-6); // ≈ 328.843
		expect(() => calculateTm(DH, DS, 1, 1, false, sCloseToZero)).toThrow(
			/Denominator is ~0|infinite Tm/i
		);
	});

	test('non-equimolar case responds to saltCorrection (A=2 µM, B=0.5 µM)', () => {
		const tm = calculateTm(DH, DS, 2, 0.5, false, 10);
		expect(tm).toBeCloseTo(42.96, 2);
	});
});

describe('isSelfComplimentary()', () => {
	// ──────────────────────────────────────────────────────────────────
	// Parameter validation
	// ──────────────────────────────────────────────────────────────────
	const badInputs = [
		['null', null],
		['undefined', undefined],
		['empty', ''],
		['number', 123],
		['array', ['A', 'T']],
		['object', { s: 'AT' }],
		['lowercase', 'atcg'],
		['mixed case', 'AtCG'],
		['invalid char', 'ATN'],
	];

	for (const [label, v] of badInputs) {
		test(`throws for invalid DNA sequence (${label})`, () => {
			// @ts-expect-error – intentionally wrong types for testing
			expect(() => isSelfComplimentary(v)).toThrow(
				/Invalid DNA sequence/i
			);
		});
	}

	// ──────────────────────────────────────────────────────────────────
	// True cases (even-length palindromic / self-complementary)
	// ──────────────────────────────────────────────────────────────────
	test('returns true for classic palindromic restriction sites', () => {
		// EcoRI
		expect(isSelfComplimentary('GAATTC')).toBe(true);
		// BamHI
		expect(isSelfComplimentary('GGATCC')).toBe(true);
		// SmaI
		expect(isSelfComplimentary('CCCGGG')).toBe(true);
		// Simple even palindrome
		expect(isSelfComplimentary('ACGT')).toBe(true);
		// Two-base pairs that are complements reversed
		expect(isSelfComplimentary('AT')).toBe(true);
		expect(isSelfComplimentary('TA')).toBe(true);
		expect(isSelfComplimentary('GC')).toBe(true);
		expect(isSelfComplimentary('CG')).toBe(true);
	});

	// ──────────────────────────────────────────────────────────────────
	// False cases
	// ──────────────────────────────────────────────────────────────────
	test('returns false for non-palindromic sequences', () => {
		expect(isSelfComplimentary('GATTACA')).toBe(false);
		expect(isSelfComplimentary('AGCTTC')).toBe(false);
		expect(isSelfComplimentary('AAAAAA')).toBe(false);
		expect(isSelfComplimentary('ATGCGA')).toBe(false);
	});

	test('returns false for single-nucleotide sequences', () => {
		expect(isSelfComplimentary('A')).toBe(false);
		expect(isSelfComplimentary('C')).toBe(false);
		expect(isSelfComplimentary('G')).toBe(false);
		expect(isSelfComplimentary('T')).toBe(false);
	});

	test('returns false for all odd-length sequences (no base self-complements)', () => {
		expect(isSelfComplimentary('ATG')).toBe(false);
		expect(isSelfComplimentary('ATGTA')).toBe(false);
		expect(isSelfComplimentary('GAATTCG')).toBe(false);
	});

	// ──────────────────────────────────────────────────────────────────
	// Consistency with reverseComplement (sanity property)
	// ──────────────────────────────────────────────────────────────────
	test('agrees with equality to reverseComplement(seq)', () => {
		const samples = [
			'GAATTC', // true
			'GGATCC', // true
			'CCCGGG', // true
			'ACGT', // true
			'GATTACA', // false
			'ATGC', // false
		];
		for (const s of samples) {
			const rc = reverseComplement(s);
			expect(isSelfComplimentary(s)).toBe(s === rc);
		}
	});

	// ──────────────────────────────────────────────────────────────────
	// Constructed palindromes: generate and test
	// ──────────────────────────────────────────────────────────────────
	test('constructed even-length palindromes always return true', () => {
		// Build some palindromes of varying even lengths by mirroring a random half
		const halves = ['AT', 'GC', 'AGTC', 'GCGC', 'AAGGTTCC'];
		for (const half of halves) {
			// Mirror via reverse-complement of the half
			const pal = half + reverseComplement(half);
			// Sanity: pal should equal its own reverse complement
			expect(reverseComplement(pal)).toBe(pal);
			// Function under test
			expect(isSelfComplimentary(pal)).toBe(true);
		}
	});

	// ──────────────────────────────────────────────────────────────────
	// Boundary: minimal even-length that is NOT self-complementary
	// ──────────────────────────────────────────────────────────────────
	test('two bases that are not reverse-complement pairs return false', () => {
		// e.g., "AA" reverse-complement is "TT" → not equal
		expect(isSelfComplimentary('AA')).toBe(false);
		expect(isSelfComplimentary('CC')).toBe(false);
		expect(isSelfComplimentary('GG')).toBe(false);
		expect(isSelfComplimentary('TT')).toBe(false);
	});
});

describe('parseThermoParamsFromResponse()', () => {
	const HTML_OK = `
	  <html><body>
	    <dH> -148000.0 </dH>
	    <dS> -410.0 </dS>
	    <saltCorrection> -11.32341669897858 </saltCorrection>
	  </body></html>`;

	test('parses all fields and converts dH to kcal/mol', () => {
		const out = parseThermoParamsFromResponse(HTML_OK);
		expect(out.dH).toBeCloseTo(-148.0, 6); // kcal/mol
		expect(out.dS).toBeCloseTo(-410.0, 6);
		expect(out.saltCorrection).toBeCloseTo(-11.32341669897858, 12);
	});

	test('tolerates whitespace/newlines around numbers', () => {
		const html = `
		  <html><body>
		    <dH>
		      -100000
		    </dH>
		    <dS>
		      -300.0
		    </dS>
		    <saltCorrection>
		      1.2345
		    </saltCorrection>
		  </body></html>`;
		const out = parseThermoParamsFromResponse(html);
		expect(out.dH).toBeCloseTo(-100.0, 6);
		expect(out.dS).toBeCloseTo(-300.0, 6);
		expect(out.saltCorrection).toBeCloseTo(1.2345, 6);
	});

	test('throws if rawHtml is not a string', () => {
		// @ts-expect-error intentional
		expect(() => parseThermoParamsFromResponse(null)).toThrow(
			/rawHtml must be a string/i
		);
	});

	test('throws if rawHtml is empty', () => {
		expect(() => parseThermoParamsFromResponse('')).toThrow(
			/non-empty string/i
		);
	});

	test('throws if dH tag missing', () => {
		const html = `<html><body><dS>-1</dS><saltCorrection>0</saltCorrection></body></html>`;
		expect(() => parseThermoParamsFromResponse(html)).toThrow(
			/dH not found|unparsable/i
		);
	});

	test('throws if dS is not numeric', () => {
		const html = `<html><body><dH>-1000</dH><dS>abc</dS><saltCorrection>1</saltCorrection></body></html>`;
		expect(() => parseThermoParamsFromResponse(html)).toThrow(
			/dS not found|unparsable/i
		);
	});

	test('throws if saltCorrection is missing', () => {
		const html = `<html><body><dH>-1000</dH><dS>-10</dS></body></html>`;
		expect(() => parseThermoParamsFromResponse(html)).toThrow(
			/saltCorrection not found|unparsable/i
		);
	});
});

describe('getThermoParams() — parameter checking', () => {
	// ────────────────────────────────────────────────────────────────
	// Sequence validation
	// ────────────────────────────────────────────────────────────────
	test('throws if seq is empty', async () => {
		await expect(() => getThermoParams('', 1, 1)).rejects.toThrow(
			/invalid dna sequence/i
		);
	});

	test('throws if seq is lowercase', async () => {
		await expect(() => getThermoParams('atgc', 1, 1)).rejects.toThrow(
			/invalid dna sequence/i
		);
	});

	test('throws if seq has invalid characters', async () => {
		await expect(() => getThermoParams('ATGXAT', 1, 1)).rejects.toThrow(
			/invalid dna sequence/i
		);
	});

	test('throws if seq is not a string (number)', async () => {
		await expect(() => getThermoParams(12345, 1, 1)).rejects.toThrow(
			/invalid dna sequence/i
		);
	});

	test('throws if seq is not a string (array)', async () => {
		await expect(() => getThermoParams(['A', 'T'], 1, 1)).rejects.toThrow(
			/invalid dna sequence/i
		);
	});

	test('throws if seq is null', async () => {
		await expect(() => getThermoParams(null, 1, 1)).rejects.toThrow(
			/invalid dna sequence/i
		);
	});

	test('throws if seq is undefined', async () => {
		await expect(() => getThermoParams(undefined, 1, 1)).rejects.toThrow(
			/invalid dna sequence/i
		);
	});

	test('throws if seq includes whitespace', async () => {
		await expect(() => getThermoParams('ATG C', 1, 1)).rejects.toThrow(
			/invalid dna sequence/i
		);
	});

	// ────────────────────────────────────────────────────────────────
	// Concentration validation (concentration)
	// ────────────────────────────────────────────────────────────────
	for (const [label, v] of [
		['zero', 0],
		['negative', -1],
		['NaN', NaN],
		['Infinity', Infinity],
		['null', null],
		['undefined', undefined],
		['non-number', '1.0'],
	]) {
		test(`throws if concentration is ${label}`, async () => {
			await expect(() => getThermoParams('ATGCT', v, 1)).rejects.toThrow(
				/concentration.*positive.*µm/i
			);
		});
	}

	// ────────────────────────────────────────────────────────────────
	// Concentration validation (limitingConc)
	// ────────────────────────────────────────────────────────────────
	for (const [label, v] of [
		['zero', 0],
		['negative', -2],
		['NaN', NaN],
		['Infinity', Infinity],
		['null', null],
		['undefined', undefined],
		['non-number', '0.5'],
	]) {
		test(`throws if limitingConc is ${label}`, async () => {
			await expect(() => getThermoParams('ATGCT', 1, v)).rejects.toThrow(
				/limitingconc.*positive.*µm/i
			);
		});
	}

	// ────────────────────────────────────────────────────────────────
	// Mismatch validation: non-object / wrong types
	// ────────────────────────────────────────────────────────────────
	test('throws if mismatch is a string', async () => {
		await expect(() =>
			getThermoParams('ATGCTAGC', 1, 1, 'oops')
		).rejects.toThrow(/invalid mismatch object/i);
	});

	test('throws if mismatch is a number', async () => {
		await expect(() =>
			getThermoParams('ATGCTAGC', 1, 1, 42)
		).rejects.toThrow(/invalid mismatch object/i);
	});

	test('throws if mismatch is an array', async () => {
		await expect(() =>
			getThermoParams('ATGCTAGC', 1, 1, [{ position: 2, type: 'A' }])
		).rejects.toThrow(/invalid mismatch object/i);
	});

	// ────────────────────────────────────────────────────────────────
	// Mismatch validation: missing / extra fields
	// ────────────────────────────────────────────────────────────────
	test('throws if mismatch missing position', async () => {
		await expect(() =>
			getThermoParams('ATGCTAGC', 1, 1, { type: 'A' })
		).rejects.toThrow(/invalid mismatch object/i);
	});

	test('throws if mismatch missing type', async () => {
		await expect(() =>
			getThermoParams('ATGCTAGC', 1, 1, { position: 2 })
		).rejects.toThrow(/invalid mismatch object/i);
	});

	test('throws if mismatch has extra keys', async () => {
		await expect(() =>
			getThermoParams('ATGCTAGC', 1, 1, {
				position: 2,
				type: 'A',
				extra: true,
			})
		).rejects.toThrow(/invalid mismatch object/i);
	});

	// ────────────────────────────────────────────────────────────────
	// Mismatch validation: position
	// ────────────────────────────────────────────────────────────────
	test('throws if mismatch.position is negative', async () => {
		await expect(() =>
			getThermoParams('ATGCTAGC', 1, 1, { position: -1, type: 'A' })
		).rejects.toThrow(/invalid mismatch object/i);
	});

	test('throws if mismatch.position not integer', async () => {
		await expect(() =>
			getThermoParams('ATGCTAGC', 1, 1, { position: 2.5, type: 'A' })
		).rejects.toThrow(/invalid mismatch object/i);
	});

	test('throws if mismatch.position equals seq.length', async () => {
		const s = 'ATGCT'; // len 5; valid last index is 4
		await expect(() =>
			getThermoParams(s, 1, 1, { position: 5, type: 'A' })
		).rejects.toThrow(/exceeds sequence length/i);
	});

	test('throws if mismatch.position exceeds seq.length', async () => {
		const s = 'ATGCT'; // len 5
		await expect(() =>
			getThermoParams(s, 1, 1, { position: 99, type: 'A' })
		).rejects.toThrow(/exceeds sequence length/i);
	});

	// ────────────────────────────────────────────────────────────────
	// Mismatch validation: type (base)
	// ────────────────────────────────────────────────────────────────
	test('throws if mismatch.type is invalid base', async () => {
		await expect(() =>
			getThermoParams('ATGCTAGC', 1, 1, { position: 2, type: 'Z' })
		).rejects.toThrow(/invalid mismatch object/i);
	});

	test('throws if mismatch.type is lowercase', async () => {
		await expect(() =>
			getThermoParams('ATGCTAGC', 1, 1, { position: 2, type: 'g' })
		).rejects.toThrow(/invalid mismatch object/i);
	});

	test('throws if mismatch.type length > 1', async () => {
		await expect(() =>
			getThermoParams('ATGCTAGC', 1, 1, { position: 2, type: 'AG' })
		).rejects.toThrow(/invalid mismatch object/i);
	});

	// ────────────────────────────────────────────────────────────────
	// Combination cases (first failing guard should trigger)
	// ────────────────────────────────────────────────────────────────
	test('sequence check occurs before concentration checks', async () => {
		await expect(() => getThermoParams('bad$x', 0, 0)).rejects.toThrow(
			/invalid dna sequence/i
		);
	});

	test('concentration check occurs before mismatch checks', async () => {
		await expect(() =>
			getThermoParams('ATGCTAGC', 0, 1, { position: 2, type: 'A' })
		).rejects.toThrow(/concentration.*positive.*µm/i);
	});

	test('limitingConc check occurs before mismatch checks', async () => {
		await expect(() =>
			getThermoParams('ATGCTAGC', 1, 0, { position: 2, type: 'A' })
		).rejects.toThrow(/limitingconc.*positive.*µm/i);
	});

	// ────────────────────────────────────────────────────────────────
	// Safe-noop mismatch inputs that should still be validated later
	// (We do NOT assert non-throw on valid inputs to avoid fetch().)
	// ────────────────────────────────────────────────────────────────
	test('null mismatch is accepted as "no mismatch" (later checks may fail)', async () => {
		await expect(() =>
			getThermoParams('ATGCTAGC', 0, 1, null)
		).rejects.toThrow(/concentration.*positive.*µm/i);
	});

	test('undefined mismatch is treated as "no mismatch" (later checks may fail)', async () => {
		await expect(() =>
			getThermoParams('ATGCTAGC', 1, 0, undefined)
		).rejects.toThrow(/limitingconc.*positive.*µm/i);
	});
});
