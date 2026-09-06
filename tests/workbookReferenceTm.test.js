import {
	calculateSnapbackTmSantaLucia,
	calculateSnapbackTmWittwer,
	calculateSnapbackTmWittwerFromStructure,
	calculateWittwerTmFromStemTm,
} from '../src/js/tm/snapbackTm.js';
import {
	HUGH_REFIT_ROWS,
	LEGACY_ALL_PEAK_SUMMARY,
	LEGACY_DELTA_TM_SUMMARY,
	LEGACY_WORKBOOK_ROWS,
	LEGACY_WORKBOOK_SOURCE,
	UNRECONSTRUCTABLE_END_ROWS,
} from './fixtures/legacySnapbackWorkbook.js';

const COMPLEMENT = Object.freeze({ A: 'T', C: 'G', G: 'C', T: 'A' });
const HISTORICAL_MONOVALENT_MM = 20;
const HISTORICAL_DNTP_MM = 0.8;
const HISTORICAL_OLIGO_UM = 0.5;
const HISTORICAL_COMPLEMENT_UM = 0.125;

const WORKBOOK_METHODS = Object.freeze([
	'empiricalNoEnds',
	'empiricalWithEnds',
	'santaLucia',
	'rochester',
]);

function reverseComplement(sequence) {
	return [...sequence]
		.reverse()
		.map((base) => COMPLEMENT[base])
		.join('');
}

function calculateHistoricalFreeMagnesiumMm(totalMagnesiumMm) {
	const totalMagnesiumM = totalMagnesiumMm / 1000;
	const totalDntpM = HISTORICAL_DNTP_MM / 1000;
	const associationConstant = 30000;
	const a = associationConstant;
	const b =
		1 + associationConstant * (totalDntpM - totalMagnesiumM);
	const c = -totalMagnesiumM;
	return (
		((-b + Math.sqrt(b ** 2 - 4 * a * c)) / (2 * a)) * 1000
	);
}

function parseWorkbookCore(coreSequence) {
	const match = coreSequence.match(
		/^([ACGT]+)\(([ACGT])([ACGT])\)([ACGT]+)$/,
	);
	if (!match) throw new Error(`Invalid workbook core sequence: ${coreSequence}`);

	const [, prefix, wildBase, variantBase, suffix] = match;
	const fullWild = `${prefix}${wildBase}${suffix}`;
	const fullVariant = `${prefix}${variantBase}${suffix}`;
	return Object.freeze({
		wildStem: fullWild.slice(1, -1),
		variantStem: fullVariant.slice(1, -1),
		snvIndex: prefix.length - 1,
		wildBase,
		variantBase,
		loopTopOutside: fullWild[0],
		tailTopOutside: fullWild.at(-1),
	});
}

function historicalConditions(row) {
	return Object.freeze({
		magnesiumMm: calculateHistoricalFreeMagnesiumMm(
			row.reportedMagnesiumMm,
		),
		monovalentMm: HISTORICAL_MONOVALENT_MM,
		concentrationUm: HISTORICAL_OLIGO_UM,
		limitingConcentrationUm: HISTORICAL_COMPLEMENT_UM,
		wittwerLogBase: 'log10',
	});
}

function hasReconstructableTerminalContexts(row) {
	const core = parseWorkbookCore(row.coreSequence);
	return (
		COMPLEMENT[core.loopTopOutside] !== row.leftMismatch &&
		COMPLEMENT[core.tailTopOutside] !== row.rightMismatch
	);
}

function buildHistoricalStructure(row) {
	const core = parseWorkbookCore(row.coreSequence);
	const tailBaseAtSnv = COMPLEMENT[core.wildBase];
	return {
		fivePrimerLimSnapExtMismatches: row.rightMismatch,
		fivePrimeStem: reverseComplement(core.wildStem),
		fivePrimeInnerLoopMismatches: '',
		// The legacy workbook records N but not loop sequence. Only N and the two
		// terminal loop bases enter the selected in-house model, so neutral filler
		// bases preserve every recorded thermodynamic input.
		stuffBetween: `${row.leftMismatch}${'A'.repeat(
			row.loopLength - 2,
		)}${core.loopTopOutside}`,
		threePrimeInnerLoopMismatches: '',
		threePrimeStem: core.wildStem,
		threePrimerLimSnapExtMismatches: core.tailTopOutside,
		threePrimerRestOfAmplicon: 'A',
		snvOnThreePrimeStem: {
			indexInThreePrimeStem: core.snvIndex,
			wildBase: core.wildBase,
			variantBase: core.variantBase,
		},
		snvOnFivePrimeStem: {
			indexInFivePrimeStem: core.wildStem.length - 1 - core.snvIndex,
			tailBaseAtSNV: tailBaseAtSnv,
			matchesWild: true,
			matchesVariant: false,
			compWildBase: tailBaseAtSnv,
			compVariantBase: COMPLEMENT[core.variantBase],
		},
	};
}

function errorStats(records) {
	const errors = records.map(({ actual, expected }) => actual - expected);
	return Object.freeze({
		n: errors.length,
		mae: errors.reduce((sum, error) => sum + Math.abs(error), 0) /
			errors.length,
		maxAbsoluteError: Math.max(...errors.map(Math.abs)),
	});
}

function summarizeErrors(errors) {
	const mean = errors.reduce((sum, value) => sum + value, 0) / errors.length;
	return Object.freeze({
		n: errors.length,
		bias: mean,
		mae:
			errors.reduce((sum, value) => sum + Math.abs(value), 0) /
			errors.length,
		rmse: Math.sqrt(
			errors.reduce((sum, value) => sum + value ** 2, 0) / errors.length,
		),
		sampleSd: Math.sqrt(
			errors.reduce((sum, value) => sum + (value - mean) ** 2, 0) /
				(errors.length - 1),
		),
	});
}

function rSquared(first, second) {
	const firstMean = first.reduce((sum, value) => sum + value, 0) / first.length;
	const secondMean =
		second.reduce((sum, value) => sum + value, 0) / second.length;
	let covariance = 0;
	let firstSquares = 0;
	let secondSquares = 0;
	for (let index = 0; index < first.length; index += 1) {
		const firstCentered = first[index] - firstMean;
		const secondCentered = second[index] - secondMean;
		covariance += firstCentered * secondCentered;
		firstSquares += firstCentered ** 2;
		secondSquares += secondCentered ** 2;
	}
	return covariance ** 2 / (firstSquares * secondSquares);
}

function expectCompatibilityEnvelope(records, { mae, maxAbsoluteError }) {
	const violations = records
		.map(({ row, actual, expected }) => ({
			workbookRow: row.workbookRow,
			error: actual - expected,
		}))
		.filter(({ error }) => Math.abs(error) > maxAbsoluteError);
	expect(violations).toEqual([]);
	const stats = errorStats(records);
	expect(stats.mae).toBeLessThanOrEqual(mae);
	expect(stats.maxAbsoluteError).toBeLessThanOrEqual(maxAbsoluteError);
}

describe('frozen Carl workbook provenance', () => {
	test('preserves every source row and identifies the two incomplete end contexts', () => {
		expect(LEGACY_WORKBOOK_ROWS).toHaveLength(47);
		expect(LEGACY_WORKBOOK_ROWS[0].workbookRow).toBe(2);
		expect(LEGACY_WORKBOOK_ROWS.at(-1).workbookRow).toBe(48);
		expect(HUGH_REFIT_ROWS).toHaveLength(44);
		expect(LEGACY_WORKBOOK_SOURCE.inputSha256).toMatch(/^[a-f0-9]{64}$/);

		const unreconstructableRows = LEGACY_WORKBOOK_ROWS.filter(
			(row) => !hasReconstructableTerminalContexts(row),
		).map((row) => row.workbookRow);
		expect(unreconstructableRows).toEqual(
			Object.keys(UNRECONSTRUCTABLE_END_ROWS).map(Number),
		);
		expect(calculateHistoricalFreeMagnesiumMm(3)).toBeCloseTo(
			2.211877,
			6,
		);
	});

	test('matches the workbook delta-Tm aggregate arithmetic', () => {
		for (const method of WORKBOOK_METHODS) {
			const errors = LEGACY_WORKBOOK_ROWS.filter(
				(row) => row.measured.matched !== null,
			).map((row) => {
				const measuredDelta =
					row.measured.matched - row.measured.mismatched;
				const predictedDelta =
					row.legacy[method].matched - row.legacy[method].mismatched;
				return predictedDelta - measuredDelta;
			});
			const actual = summarizeErrors(errors);
			const expected = LEGACY_DELTA_TM_SUMMARY[method];
			expect(actual.n).toBe(expected.n);
			expect(actual.bias).toBeCloseTo(expected.bias, 12);
			expect(actual.mae).toBeCloseTo(expected.mae, 12);
			expect(actual.rmse).toBeCloseTo(expected.rmse, 12);
			expect(actual.sampleSd).toBeCloseTo(expected.sampleSd, 12);
		}
	});

	test('matches the workbook all-peak summary arithmetic', () => {
		for (const method of WORKBOOK_METHODS) {
			const measured = [];
			const predicted = [];
			for (const row of LEGACY_WORKBOOK_ROWS) {
				for (const peak of ['matched', 'mismatched']) {
					if (row.measured[peak] === null) continue;
					measured.push(row.measured[peak]);
					predicted.push(row.legacy[method][peak]);
				}
			}
			const actual = summarizeErrors(
				predicted.map((value, index) => value - measured[index]),
			);
			const expected = LEGACY_ALL_PEAK_SUMMARY[method];
			expect(actual.n).toBe(expected.n);
			expect(actual.bias).toBeCloseTo(expected.bias, 12);
			expect(actual.mae).toBeCloseTo(expected.mae, 12);
			expect(actual.sampleSd).toBeCloseTo(expected.sampleSd, 12);
			expect(rSquared(measured, predicted)).toBeCloseTo(
				expected.rSquared,
				4,
			);
		}
	});

	test('preserves the workbook response to increased magnesium', () => {
		const bySequence = new Map();
		for (const row of LEGACY_WORKBOOK_ROWS) {
			const key = [
				row.leftMismatch,
				row.coreSequence,
				row.rightMismatch,
				row.loopLength,
			].join('|');
			if (!bySequence.has(key)) bySequence.set(key, []);
			bySequence.get(key).push(row);
		}
		const titrations = [...bySequence.values()].filter(
			(rows) => new Set(rows.map((row) => row.reportedMagnesiumMm)).size > 1,
		);
		expect(titrations).toHaveLength(14);

		for (const rows of titrations) {
			const byMagnesium = new Map();
			for (const row of rows) {
				if (!byMagnesium.has(row.reportedMagnesiumMm)) {
					byMagnesium.set(row.reportedMagnesiumMm, row);
				}
			}
			const ordered = [...byMagnesium.values()].sort(
				(left, right) =>
					left.reportedMagnesiumMm - right.reportedMagnesiumMm,
			);
			for (let index = 1; index < ordered.length; index += 1) {
				for (const method of WORKBOOK_METHODS) {
					for (const peak of ['matched', 'mismatched']) {
						expect(ordered[index].legacy[method][peak]).toBeGreaterThan(
							ordered[index - 1].legacy[method][peak],
						);
					}
				}
			}
		}
	});

	test('preserves identical predictions for duplicate calculation inputs', () => {
		const groups = new Map();
		for (const row of LEGACY_WORKBOOK_ROWS) {
			const key = [
				row.leftMismatch,
				row.coreSequence,
				row.rightMismatch,
				row.reportedMagnesiumMm,
				row.loopLength,
			].join('|');
			if (!groups.has(key)) groups.set(key, []);
			groups.get(key).push(row);
		}
		const duplicates = [...groups.values()].filter((rows) => rows.length > 1);
		expect(duplicates).toHaveLength(6);
		for (const rows of duplicates) {
			for (const row of rows.slice(1)) {
				expect(row.legacy).toEqual(rows[0].legacy);
			}
		}
	});
});

describe('Carl workbook empirical composition', () => {
	test.each(HUGH_REFIT_ROWS)(
		'round-trips $snpName $peak from Hugh Refit row $workbookRow',
		({ recoveredStemTm, loopLength, originalEmpiricalTm }) => {
			const actual = calculateWittwerTmFromStemTm(
				recoveredStemTm,
				loopLength,
				'log10',
			);
			expect(actual).toBeCloseTo(originalEmpiricalTm, 10);
		},
	);

	test('reproduces empirical no-ends benchmarks across all 47 experiments', () => {
		const matched = [];
		const mismatched = [];
		for (const row of LEGACY_WORKBOOK_ROWS) {
			const core = parseWorkbookCore(row.coreSequence);
			const conditions = historicalConditions(row);
			matched.push({
				row,
				actual: calculateSnapbackTmWittwer(
					core.wildStem,
					row.loopLength,
					undefined,
					conditions,
				),
				expected: row.legacy.empiricalNoEnds.matched,
			});
			mismatched.push({
				row,
				actual: calculateSnapbackTmWittwer(
					core.variantStem,
					row.loopLength,
					{
						position: core.snvIndex,
						type: COMPLEMENT[core.wildBase],
					},
					conditions,
				),
				expected: row.legacy.empiricalNoEnds.mismatched,
			});
		}

		expectCompatibilityEnvelope(matched, {
			mae: 0.2,
			maxAbsoluteError: 0.5,
		});
		expectCompatibilityEnvelope(mismatched, {
			mae: 0.25,
			maxAbsoluteError: 0.5,
		});
	});

	test('reproduces reconstructable terminal-end empirical benchmarks', () => {
		const matched = [];
		const mismatched = [];
		const rows = LEGACY_WORKBOOK_ROWS.filter(
			hasReconstructableTerminalContexts,
		);
		expect(rows).toHaveLength(45);

		for (const row of rows) {
			const actual = calculateSnapbackTmWittwerFromStructure(
				buildHistoricalStructure(row),
				historicalConditions(row),
			);
			matched.push({
				row,
				actual: actual.wildTm,
				expected: row.legacy.empiricalWithEnds.matched,
			});
			mismatched.push({
				row,
				actual: actual.variantTm,
				expected: row.legacy.empiricalWithEnds.mismatched,
			});
		}

		expectCompatibilityEnvelope(matched, {
			mae: 0.25,
			maxAbsoluteError: 1,
		});
		expectCompatibilityEnvelope(mismatched, {
			mae: 0.35,
			maxAbsoluteError: 1,
		});
	});
});

describe('Carl workbook SantaLucia compatibility', () => {
	test('reproduces reconstructable benchmarks within the historical-model envelope', () => {
		const matched = [];
		const mismatched = [];
		const rows = LEGACY_WORKBOOK_ROWS.filter(
			hasReconstructableTerminalContexts,
		);
		expect(rows).toHaveLength(45);

		for (const row of rows) {
			const actual = calculateSnapbackTmSantaLucia(
				buildHistoricalStructure(row),
				historicalConditions(row),
			);
			matched.push({
				row,
				actual: actual.wildTm,
				expected: row.legacy.santaLucia.matched,
			});
			mismatched.push({
				row,
				actual: actual.variantTm,
				expected: row.legacy.santaLucia.mismatched,
			});
		}

		expectCompatibilityEnvelope(matched, {
			mae: 0.65,
			maxAbsoluteError: 1.5,
		});
		expectCompatibilityEnvelope(mismatched, {
			mae: 0.55,
			maxAbsoluteError: 1.5,
		});
	});
});
