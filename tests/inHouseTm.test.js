import {
	calculateDuplexThermodynamics,
	calculateDuplexThermodynamicsFromStrands,
	calculateSnapbackTmSantaLucia,
	calculateSnapbackTmWittwer,
	calculateSnapbackTmWittwerFromStructure,
	getSantaLuciaHairpinLoopParams,
} from '../src/js/tm/snapbackTm.js';
import { DNA_COMPLEMENT } from '../src/js/tm/parameters.js';
import { calculateOwczarzySaltCorrection } from '../src/js/tm/saltCorrection.js';
import { readStoredTmConditions } from '../src/js/shared/tmConditions.js';
import {
	renderDeltaTmTable,
	renderTmSummary,
} from '../src/js/pages/resultsRender.js';
import {
	calculateSnapbackTmSantaLucia as calculateSnapbackTmSantaLuciaPublic,
	calculateSnapbackTmWittwer as calculateSnapbackTmWittwerPublic,
	createSnapback,
} from '../dist/script.js';

const conditions = { magnesiumMm: 3, monovalentMm: 13.7 };

function naturalLoopStructure() {
	return {
		fivePrimerLimSnapExtMismatches: 'A',
		fivePrimeStem: 'GGCTGGC',
		fivePrimeInnerLoopMismatches: '',
		stuffBetween: 'ATACGA',
		threePrimeInnerLoopMismatches: '',
		threePrimeStem: 'GCCAGCC',
		threePrimerLimSnapExtMismatches: 'A',
		threePrimerRestOfAmplicon: 'CGT',
		snvOnThreePrimeStem: {
			indexInThreePrimeStem: 3,
			wildBase: 'A',
			variantBase: 'T',
		},
		snvOnFivePrimeStem: {
			indexInFivePrimeStem: 3,
			tailBaseAtSNV: 'T',
			matchesWild: true,
			matchesVariant: false,
			compWildBase: 'T',
			compVariantBase: 'A',
		},
	};
}

describe('in-house stem thermodynamics', () => {
	test('sums SantaLucia/Hicks initiation, stack, and both terminal-AT terms', () => {
		const result = calculateDuplexThermodynamics('AA');
		expect(result.dH).toBeCloseTo(-3.0, 12);
		expect(result.dS).toBeCloseTo(-13.2, 12);
		expect(result.terminalAT).toEqual(['left', 'right']);
	});

	test('matches the saved Carl LabVIEW matched-sequence H/S vector', () => {
		const result = calculateDuplexThermodynamics('AAGTGCACTTAGCTACGTAG');
		expect(result.dH).toBeCloseTo(-154.2, 12);
		expect(result.dS).toBeCloseTo(-421.3, 12);
	});

	test('matches the saved Carl LabVIEW internal-mismatch H/S vector', () => {
		const result = calculateDuplexThermodynamics(
			'AAGTGCACTTAGCTACGTAG',
			{ position: 3, type: 'G' },
		);
		expect(result.dH).toBeCloseTo(-145.8, 12);
		expect(result.dS).toBeCloseTo(-400.2, 12);
	});

	test('uses the two mismatch tetrads adjacent to a single internal mismatch', () => {
		const result = calculateDuplexThermodynamicsFromStrands('AGA', 'TTT');
		expect(result.mismatchPositions).toEqual([1]);
		expect(result.stacks.map((stack) => stack.key)).toEqual([
			'AG/TT',
			'GA/TT',
		]);
		expect(result.stacks.every((stack) => stack.kind === 'internal-mismatch')).toBe(
			true,
		);
	});

	test('has parameters for every orientation around every single noncanonical pair', () => {
		for (const topMismatch of 'ACGT') {
			for (const bottomMismatch of 'ACGT') {
				if (DNA_COMPLEMENT[topMismatch] === bottomMismatch) continue;
				for (const left of 'ACGT') {
					for (const right of 'ACGT') {
						const top = `${left}${topMismatch}${right}`;
						const bottomAligned = `${DNA_COMPLEMENT[left]}${bottomMismatch}${DNA_COMPLEMENT[right]}`;
						expect(() =>
							calculateDuplexThermodynamicsFromStrands(
								top,
								[...bottomAligned].reverse().join(''),
							),
						).not.toThrow();
					}
				}
			}
		}
	});
});

describe('loop and salt models', () => {
	test('does not reuse an invalid saved zero-ion pair for primer previews', () => {
		const storage = {
			getItem: (key) => ({ magnesiumMm: '0', monovalentMm: '0' })[key],
		};
		expect(readStoredTmConditions(storage)).toEqual({
			magnesiumMm: 2.2,
			monovalentMm: 13.7,
		});
	});

	test('uses the corrected 8- and 10-base loop anchors', () => {
		expect(getSantaLuciaHairpinLoopParams(8).dS).toBe(-13.9);
		expect(getSantaLuciaHairpinLoopParams(10).dS).toBe(-14.8);
	});

	test('interpolates an omitted loop size in free-energy space', () => {
		const loop = getSantaLuciaHairpinLoopParams(11);
		expect(loop.dG37).toBeCloseTo(4.8, 12);
		expect(loop.dS).toBeCloseTo((-4.8 * 1000) / 310.15, 12);
		expect(loop.interpolated).toBe(true);
	});

	test('reproduces the Owczarzy paper correction factor example', () => {
		const result = calculateOwczarzySaltCorrection({
			stemSequence: 'AAGGCGAGTCAGGCTCAGTG',
			stemDeltaH: -160,
			magnesiumMm: 1.5,
			monovalentMm: 5,
		});
		expect(result.ratio).toBeCloseTo(7.74597, 5);
		expect(result.correctionFactor).toBeCloseTo(7.0537e-5, 9);
		expect(result.regime).toBe('magnesium-dominant');
	});

	test('covers every Owczarzy ion regime and rejects the undefined zero-ion case', () => {
		const base = {
			stemSequence: 'AAGGCGAGTCAGGCTCAGTG',
			stemDeltaH: -160,
		};
		expect(
			calculateOwczarzySaltCorrection({
				...base,
				magnesiumMm: 0,
				monovalentMm: 50,
			}).regime,
		).toBe('monovalent-only');
		expect(
			calculateOwczarzySaltCorrection({
				...base,
				magnesiumMm: 0.01,
				monovalentMm: 100,
			}).regime,
		).toBe('monovalent-dominant');
		expect(
			calculateOwczarzySaltCorrection({
				...base,
				magnesiumMm: 1,
				monovalentMm: 50,
			}).regime,
		).toBe('mixed-ion');
		expect(
			calculateOwczarzySaltCorrection({
				...base,
				magnesiumMm: 1,
				monovalentMm: 0,
			}).regime,
		).toBe('magnesium-only');
		expect(() =>
			calculateOwczarzySaltCorrection({
				...base,
				magnesiumMm: 0,
				monovalentMm: 0,
			}),
		).toThrow(/requires magnesiumMm or monovalentMm/);
	});
});

describe('complete snapback methods', () => {
	test('uses log10 exclusively for the Carl/Wittwer loop term', async () => {
		const defaultResult = calculateSnapbackTmWittwer(
			'ACGTACGT',
			6,
			undefined,
			conditions,
		);
		const log10 = calculateSnapbackTmWittwer(
			'ACGTACGT',
			6,
			undefined,
			conditions,
			{ logBase: 'log10' },
		);
		const log10ViaConditions = calculateSnapbackTmWittwer(
			'ACGTACGT',
			6,
			undefined,
			{ ...conditions, wittwerLogBase: 'log10' },
		);
		expect(defaultResult).toBe(log10);
		expect(log10ViaConditions).toBe(log10);
		expect(() =>
			calculateSnapbackTmWittwer('ACGTACGT', 6, undefined, {
				...conditions,
				wittwerLogBase: 'ln',
			}),
		).toThrow(/must be "log10"/);
		expect(() =>
			calculateSnapbackTmWittwer(
				'ACGTACGT',
				6,
				undefined,
				conditions,
				{ logBase: 'ln' },
			),
		).toThrow(/uses log10 only/);
		await expect(
			calculateSnapbackTmWittwerPublic(
				'ACGTACGT',
				6,
				undefined,
				conditions,
			),
		).resolves.toBe(defaultResult);
		await expect(
			calculateSnapbackTmWittwerPublic(
				'ACGTACGT',
				6,
				undefined,
				{ ...conditions, wittwerLogBase: 'ln' },
			),
		).rejects.toThrow(/must be "log10"/);
	});

	test('uses one natural loop mismatch and one extension-blocking mismatch', () => {
		const result = calculateSnapbackTmSantaLucia(
			naturalLoopStructure(),
			conditions,
		);
		expect(result.components.loop.N).toBe(6);
		expect(result.components.terminalMismatches).toHaveLength(2);
		expect(result.components.terminalMismatches[0].label).toBe(
			'natural-loop-side',
		);
		expect(result.components.terminalMismatches).toEqual([
			{ label: 'natural-loop-side', top2: 'AG', bottom2: 'AC' },
			{ label: 'extension-blocking-end', top2: 'CA', bottom2: 'GA' },
		]);
		expect(result.alleles.wild.terminalMismatches.entries).toEqual([
			{
				label: 'natural-loop-side',
				top2: 'AG',
				bottom2: 'AC',
				dH: -4.6,
				dS: -11.6,
			},
			{
				label: 'extension-blocking-end',
				top2: 'CA',
				bottom2: 'GA',
				dH: -4.6,
				dS: -11.6,
			},
		]);
		expect(Number.isFinite(result.wildTm)).toBe(true);
		expect(Number.isFinite(result.variantTm)).toBe(true);
	});

	test('counts an engineered loop mismatch on both sides of the loop', () => {
		const extended = naturalLoopStructure();
		extended.stuffBetween = 'ATACG';
		extended.threePrimeInnerLoopMismatches = 'T';
		extended.fivePrimeInnerLoopMismatches = 'G';
		const result = calculateSnapbackTmSantaLucia(extended, conditions);
		expect(result.components.loop.N).toBe(7);
		expect(result.components.terminalMismatches[0].label).toBe(
			'engineered-loop-side',
		);
	});

	test('requires three matched stem bases on each side of the internal SNV', () => {
		const extended = naturalLoopStructure();
		extended.snvOnThreePrimeStem.indexInThreePrimeStem = 2;
		extended.snvOnFivePrimeStem.indexInFivePrimeStem = 4;
		expect(() => calculateSnapbackTmSantaLucia(extended, conditions)).toThrow(
			/at least 3 matched stem bases on each side/,
		);
	});

	test('rejects more than one engineered mismatch at either snapback end', () => {
		const extended = naturalLoopStructure();
		extended.fivePrimerLimSnapExtMismatches = 'AA';
		expect(() => calculateSnapbackTmSantaLucia(extended, conditions)).toThrow(
			/exactly one extension-blocking mismatch base/,
		);

		const loopExtended = naturalLoopStructure();
		loopExtended.fivePrimeInnerLoopMismatches = 'GG';
		loopExtended.threePrimeInnerLoopMismatches = 'TT';
		expect(() => calculateSnapbackTmSantaLucia(loopExtended, conditions)).toThrow(
			/either one engineered mismatch pair or one natural mismatch/,
		);
	});

	test('rejects contradictory tail allele metadata', () => {
		const extended = naturalLoopStructure();
		extended.snvOnFivePrimeStem.matchesWild = false;
		expect(() => calculateSnapbackTmSantaLucia(extended, conditions)).toThrow(
			/internally inconsistent/,
		);
	});

	test('returns a complete empirical/Wittwer comparison from the same structure', () => {
		const result = calculateSnapbackTmWittwerFromStructure(
			naturalLoopStructure(),
			conditions,
		);
		expect(result.loopLength).toBe(6);
		expect(result.logBase).toBe('log10');
		expect(result.terminalMismatches).toHaveLength(2);
		expect(Number.isFinite(result.wildTm)).toBe(true);
		expect(Number.isFinite(result.variantTm)).toBe(true);
	});

	test('allows empirical strand concentrations to be changed through conditions', () => {
		const equal = calculateSnapbackTmWittwerFromStructure(
			naturalLoopStructure(),
			{ ...conditions, concentrationUm: 0.5, limitingConcentrationUm: 0.5 },
		);
		const asymmetric = calculateSnapbackTmWittwerFromStructure(
			naturalLoopStructure(),
			{ ...conditions, concentrationUm: 0.5, limitingConcentrationUm: 0.1 },
		);
		expect(asymmetric.wildTm).not.toBe(equal.wildTm);
		expect(
			asymmetric.alleles.wild.duplex.conditions.limitingConcentrationUm,
		).toBe(0.1);
	});

	test('keeps the intramolecular SantaLucia result concentration-independent', () => {
		const first = calculateSnapbackTmSantaLucia(naturalLoopStructure(), {
			...conditions,
			concentrationUm: 0.5,
			limitingConcentrationUm: 0.5,
		});
		const second = calculateSnapbackTmSantaLucia(naturalLoopStructure(), {
			...conditions,
			concentrationUm: 5,
			limitingConcentrationUm: 0.1,
		});
		expect(second.wildTm).toBe(first.wildTm);
		expect(second.variantTm).toBe(first.variantTm);
	});
});

describe('results presentation', () => {
	test('shows SantaLucia as primary and Wittwer as the comparison', () => {
		document.body.innerHTML = `
			<span id="wildTm"></span>
			<span id="varTm"></span>
			<span id="wittwerWildTm"></span>
			<span id="wittwerVarTm"></span>
		`;
		renderTmSummary({
			snapbackMeltingTms: { wildTm: 60.12, variantTm: 51.34 },
			snapbackTmWittwer: { wildTm: 58.76, variantTm: 49.87 },
		});
		expect(document.getElementById('wildTm').textContent).toBe('60.1');
		expect(document.getElementById('varTm').textContent).toBe('51.3');
		expect(document.getElementById('wittwerWildTm').textContent).toBe('58.8');
		expect(document.getElementById('wittwerVarTm').textContent).toBe('49.9');
	});
});

describe('web-app integration', () => {
	test.each([39, 39.5, 81])(
		'rejects target Tm %p through the public API just as the form does',
		async (targetTm) => {
			await expect(
				createSnapback(
					'A'.repeat(20) +
						'CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA' +
						'T'.repeat(20),
					20,
					20,
					{ index: 40, variantBase: 'A' },
					targetTm,
					{ magnesiumMm: 2.2, monovalentMm: 13.7 },
				),
			).rejects.toThrow('whole number from 40 to 80');
		},
	);

	test('does not reject a valid warm stem merely because a sub-40 stem is closer', async () => {
		const sequence =
			'ATGCGCTAACTCAGGGAGCCTGTTGCGACGTTGGGGTCCAAGTTTTATATGTTACCCTTGGCCTCCCCTAAGCCGAGAGGTTAGCTAAGCTGTCGGCGCACACCCCATATAAATGCTCCT';
		const result = await createSnapback(
			sequence,
			20,
			20,
			{ index: 60, variantBase: 'T' },
			40,
			{ magnesiumMm: 2.2, monovalentMm: 13.7 },
		);
		expect(result.snapbackMeltingTms.wildTm).toBeGreaterThanOrEqual(40);
		expect(
			result.optimizedSnapbackOptions.onReversePrimer.matchVariant.wildTm,
		).toBeCloseTo(48.16, 2);
		const selectedSide = result.tailOnForwardPrimer
			? 'onForwardPrimer'
			: 'onReversePrimer';
		const selectedMatch = result.matchesWild ? 'matchWild' : 'matchVariant';
		const selectedDelta = result.meltingTempDiffs[selectedSide][selectedMatch];
		const allDeltas = [
			...Object.values(result.meltingTempDiffs.onForwardPrimer),
			...Object.values(result.meltingTempDiffs.onReversePrimer),
		].filter(Number.isFinite);
		expect(selectedDelta).toBe(Math.max(...allDeltas));
	});

	test('keeps viable options when one cannot reach 40 C and renders that option as unavailable', async () => {
		const sequence =
			'AATATCTATGTATTCATATGGTTAGCTCTTTTTTACTATAAATATTATTTTGTCGATTAATTCTTTCAAT';
		const result = await createSnapback(
			sequence,
			28,
			28,
			{ index: 35, variantBase: 'G' },
			60,
			{ magnesiumMm: 2.2, monovalentMm: 13.7 },
		);
		expect(result.optimizedSnapbackOptions.onForwardPrimer.matchVariant).toBeNull();
		expect(result.meltingTempDiffs.onForwardPrimer.matchVariant).toBeNull();
		expect(result.tailOnForwardPrimer).toBe(false);
		expect(result.matchesWild).toBe(true);
		expect(result.snapbackTmWittwer.logBase).toBe('log10');

		const santaPublic = await calculateSnapbackTmSantaLuciaPublic(
			result.descriptiveExtendedSnapback,
			result.tmConditions,
		);
		const wittwerPublic = await calculateSnapbackTmWittwerPublic(
			result.descriptiveExtendedSnapback,
			result.tmConditions,
		);
		expect(santaPublic).toEqual(result.snapbackTmSantaLucia);
		expect(wittwerPublic).toEqual(result.snapbackTmWittwer);

		document.body.innerHTML = `
			<span id="dt-wild-heading"></span>
			<span id="dt-var-heading"></span>
			<span id="dt-fwd-heading"></span>
			<span id="dt-rev-heading"></span>
			<span id="dt-fwd-wild"></span>
			<span id="dt-fwd-var"></span>
			<span id="dt-rev-wild"></span>
			<span id="dt-rev-var"></span>
		`;
		renderDeltaTmTable(result, sequence[35], 'G');
		expect(document.getElementById('dt-fwd-var').textContent).toBe('—');
	});
});
