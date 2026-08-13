import { jest } from '@jest/globals';

import {
	buildTmRequestParams,
	createSnapback,
	getOligoTm,
} from '../dist/script.js';
import {
	renderDeltaTmTable,
	renderSnapbackPrimer,
	renderTmConditions,
} from '../src/js/pages/resultsRender.js';
import { renderStemDiagram } from '../src/js/pages/resultsStemDiagram.js';
import { validateTmConditions } from '../src/js/shared/validators.js';

describe('Carl-requested result annotations', () => {
	afterEach(() => {
		document.body.innerHTML = '';
		jest.useRealTimers();
	});

	test('shows the first naturally mismatched bases inside the loop', () => {
		document.body.innerHTML = `
			<div id="stemDiagramWrapper">
				<div id="stemDiagram"></div>
				<div id="stemSnvLabel"></div>
			</div>
		`;

		renderStemDiagram(
			{
				fivePrimeStem: 'ACGT',
				fivePrimeInnerLoopMismatches: '',
				fivePrimerLimSnapExtMismatches: '',
				forwardPrimer: 'ACGTACGTACGT',
			},
			{
				threePrimeStem: 'ACGT',
				threePrimeInnerLoopMismatches: '',
				threePrimerLimSnapExtMismatches: '',
				threePrimerRestOfAmplicon: '',
				stuffBetween: 'ACCCCCG',
			},
			'A',
			'G',
			null,
			true,
		);

		const natural = document.querySelectorAll(
			'.stem-mismatch-block--inner-loop.stem-mismatch-block--natural',
		);
		expect(natural).toHaveLength(2);
		expect(natural[0].textContent).toBe('A');
		expect(natural[1].textContent).toBe('G');
	});

	test('marks both engineered mismatch regions in the snapback primer sequence', () => {
		document.body.innerHTML = `
			<pre id="snapSeq"></pre>
			<span id="snapPrimerLabel"></span>
		`;

		renderSnapbackPrimer(
			{
				tailOnForwardPrimer: true,
				snapbackSeq: 'GACGTCTTAA',
				descriptiveUnExtendedSnapbackPrimer: {
					fivePrimerLimSnapExtMismatches: 'G',
					fivePrimeStem: 'ACGT',
					fivePrimeInnerLoopMismatches: 'C',
					forwardPrimer: 'TTAA',
				},
			},
			4,
			4,
		);

		expect(document.getElementById('snapSeq').textContent).toBe(
			'GACGTCTTAA',
		);
		expect(
			document.querySelector('.seq-seg--terminal-mismatch').textContent,
		).toBe('G');
		expect(
			document.querySelector('.seq-seg--terminal-mismatch').classList,
		).toContain('seq-seg--tail');
		expect(
			document.querySelector('.seq-seg--inner-loop-mismatch').textContent,
		).toBe('C');
		expect(
			document.querySelector('.seq-seg--inner-loop-mismatch').classList,
		).toContain('seq-seg--tail');
		expect(document.querySelector('.seq-seg--stem').textContent).toBe(
			'ACGT',
		);
		expect(document.querySelector('.seq-seg--primer').textContent).toBe(
			'TTAA',
		);
	});

	test('boxes a natural first-loop primer base in blue when no base was added', () => {
		document.body.innerHTML = `
			<pre id="snapSeq"></pre>
			<span id="snapPrimerLabel"></span>
		`;

		renderSnapbackPrimer(
			{
				tailOnForwardPrimer: true,
				snapbackSeq: 'GACGTATTAA',
				descriptiveUnExtendedSnapbackPrimer: {
					fivePrimerLimSnapExtMismatches: 'G',
					fivePrimeStem: 'ACGT',
					fivePrimeInnerLoopMismatches: '',
					forwardPrimer: 'ATTAA',
				},
				descriptiveExtendedSnapback: {
					stuffBetween: 'CCCG',
				},
			},
			5,
			5,
		);

		const natural = document.querySelector(
			'.seq-seg--inner-loop-mismatch.seq-seg--natural-mismatch',
		);
		expect(natural.textContent).toBe('A');
		expect(natural.classList).toContain('seq-seg--primer');
		expect(natural.classList).not.toContain('seq-seg--tail');
		expect(document.getElementById('snapSeq').textContent).toBe(
			'GACGTATTAA',
		);
	});

	test('keeps reverse-orientation wild and variant labels aligned', () => {
		jest.useFakeTimers();
		document.body.innerHTML = `
			<div id="stemDiagramWrapper">
				<div id="stemDiagram"></div>
				<div id="stemSnvLabel"></div>
			</div>
		`;

		renderStemDiagram(
			{
				fivePrimeStem: 'TTTTCTTT',
				fivePrimeInnerLoopMismatches: '',
				fivePrimerLimSnapExtMismatches: '',
				forwardPrimer: 'ACGTACGTACGT',
			},
			{
				threePrimeStem: 'AAAAGAAA',
				threePrimeInnerLoopMismatches: '',
				threePrimerLimSnapExtMismatches: '',
				threePrimerRestOfAmplicon: '',
				stuffBetween: 'ACCCCCG',
			},
			'T',
			'C',
			4,
			false,
		);

		expect(document.querySelector('.stem-nt--snv').textContent).toBe('G');
		expect(document.getElementById('stemSnvLabel').textContent).toBe(
			'Variant: G',
		);

		jest.advanceTimersByTime(2000);
		expect(document.querySelector('.stem-nt--snv').textContent).toBe('A');
		expect(document.getElementById('stemSnvLabel').textContent).toBe(
			'Wild: A',
		);
		jest.clearAllTimers();
	});

	test('shows the actual mismatch pair in each delta-Tm cell', () => {
		document.body.innerHTML = `
			<table>
				<thead><tr>
					<th></th>
					<th id="dt-wild-heading"></th>
					<th id="dt-var-heading"></th>
				</tr></thead>
				<tbody>
					<tr><th id="dt-fwd-heading"></th><td id="dt-fwd-wild"></td><td id="dt-fwd-var"></td></tr>
					<tr><th id="dt-rev-heading"></th><td id="dt-rev-wild"></td><td id="dt-rev-var"></td></tr>
				</tbody>
			</table>
		`;

		renderDeltaTmTable(
			{
				meltingTempDiffs: {
					onForwardPrimer: { matchWild: 10, matchVariant: 9 },
					onReversePrimer: { matchWild: 8, matchVariant: 7 },
				},
			},
			'A',
			'G',
		);

		expect(document.getElementById('dt-wild-heading').textContent).toBe(
			'Wild-type match',
		);
		expect(document.getElementById('dt-var-heading').textContent).toBe(
			'Variant match',
		);
		expect(document.getElementById('dt-fwd-heading').textContent).toBe(
			'Tail on forward primer',
		);
		expect(document.getElementById('dt-rev-heading').textContent).toBe(
			'Tail on reverse primer',
		);
		expect(document.getElementById('dt-fwd-wild').textContent).toBe(
			'10.0 (G-T)',
		);
		expect(document.getElementById('dt-fwd-var').textContent).toBe(
			'9.0 (A-C)',
		);
		expect(document.getElementById('dt-rev-wild').textContent).toBe(
			'8.0 (C-A)',
		);
		expect(document.getElementById('dt-rev-var').textContent).toBe(
			'7.0 (T-G)',
		);
	});

	test('colors natural loop bases as primer and engineered loop/end bases as tail', () => {
		document.body.innerHTML = `
			<div id="stemDiagramWrapper">
				<div id="stemDiagram"></div>
				<div id="stemSnvLabel"></div>
			</div>
		`;

		renderStemDiagram(
			{
				fivePrimeStem: 'ACGT',
				fivePrimeInnerLoopMismatches: '',
				fivePrimerLimSnapExtMismatches: 'G',
				forwardPrimer: 'ACGTACGTACGT',
			},
			{
				threePrimeStem: 'ACGT',
				threePrimeInnerLoopMismatches: '',
				threePrimerLimSnapExtMismatches: 'C',
				threePrimerRestOfAmplicon: '',
				stuffBetween: 'ACCCCCG',
			},
			'A',
			'G',
			null,
			true,
		);

		expect(
			document.querySelector('.stem-row--top .stem-mismatch-block--natural')
				.classList,
		).toContain('stem-mismatch-block--primer');
		expect(
			document.querySelector('.stem-row--top .stem-mismatch-block--terminal')
				.classList,
		).toContain('stem-mismatch-block--tail');

		renderStemDiagram(
			{
				fivePrimeStem: 'ACGT',
				fivePrimeInnerLoopMismatches: 'T',
				fivePrimerLimSnapExtMismatches: 'G',
				forwardPrimer: 'ACGTACGTACGT',
			},
			{
				threePrimeStem: 'ACGT',
				threePrimeInnerLoopMismatches: 'A',
				threePrimerLimSnapExtMismatches: 'C',
				threePrimerRestOfAmplicon: '',
				stuffBetween: 'ACCCCCG',
			},
			'A',
			'G',
			null,
			true,
		);

		expect(
			document.querySelector('.stem-row--top .stem-mismatch-block--inner-loop')
				.classList,
		).toContain('stem-mismatch-block--tail');
	});
});

describe('Carl-requested cation inputs', () => {
	afterEach(() => {
		jest.restoreAllMocks();
	});
	test('accepts the requested default free-magnesium and monovalent values', () => {
		expect(validateTmConditions('3.0', '13.7')).toEqual({
			ok: true,
			msg: '',
			data: { magnesiumMm: 3, monovalentMm: 13.7 },
		});
	});

	test('labels the exact ionic quantities used on the results page', () => {
		document.body.innerHTML = `
			<span id="freeMagnesiumMm"></span>
			<span id="totalMonovalentMm"></span>
		`;

		renderTmConditions({ magnesiumMm: 4.6, monovalentMm: 17.2 });

		expect(document.getElementById('freeMagnesiumMm').textContent).toBe(
			'4.6',
		);
		expect(document.getElementById('totalMonovalentMm').textContent).toBe(
			'17.2',
		);
	});

	test.each([
		['', '13.7'],
		['3.0', ''],
		['-0.1', '13.7'],
		['3.0', '-0.1'],
		['not-a-number', '13.7'],
	])('rejects invalid cation values (%s, %s)', (magnesium, monovalent) => {
		expect(validateTmConditions(magnesium, monovalent).ok).toBe(false);
	});

	test('builds legacy API parameters from explicit cation values', () => {
		const params = buildTmRequestParams('ACGTACGT', null, {
			magnesiumMm: 4.2,
			monovalentMm: 11.5,
		});

		expect(params.get('mg')).toBe('4.2');
		expect(params.get('mono')).toBe('11.5');
	});

	test('uses 3.0 mM magnesium and 13.7 mM monovalent by default', () => {
		const params = buildTmRequestParams('ACGTACGT');

		expect(params.get('mg')).toBe('3');
		expect(params.get('mono')).toBe('13.7');
	});

	test('passes cation values to the live-request URL', async () => {
		const fetchSpy = jest.spyOn(global, 'fetch').mockResolvedValue({
			ok: true,
			text: async () => '<html><body><tm>42.5</tm></body></html>',
		});

		await expect(
			getOligoTm('ACGTACGT', null, {
				magnesiumMm: 5.1,
				monovalentMm: 9.4,
			}),
		).resolves.toBe(42.5);

		const requested = new URL(fetchSpy.mock.calls[0][0]);
		const directUrl = requested.searchParams.get('url') || requested.href;
		const direct = new URL(directUrl);
		expect(direct.searchParams.get('mg')).toBe('5.1');
		expect(direct.searchParams.get('mono')).toBe('9.4');

		fetchSpy.mockRestore();
	});

	test('passes custom cations to every request in a complete design', async () => {
		const fetchSpy = jest.spyOn(global, 'fetch').mockImplementation((url) => {
			const outer = new URL(url);
			const direct = new URL(
				outer.searchParams.get('url') || outer.href,
			);
			const seq = direct.searchParams.get('seq') || '';
			const tm = 20 + seq.length * 2;
			const mismatchTm = tm - 10;
			return Promise.resolve({
				ok: true,
				text: async () =>
					`<html><body><tm>${tm}</tm><mmtm>${mismatchTm}</mmtm>` +
					'<dH>-80000</dH><dS>-220</dS>' +
					'<saltCorrection>-5</saltCorrection></body></html>',
			});
		});

		const sequence = 'ACGT'.repeat(25);
		await createSnapback(
			sequence,
			20,
			20,
			{ index: 50, variantBase: 'A' },
			60,
			{ magnesiumMm: 4.6, monovalentMm: 17.2 },
		);

		expect(fetchSpy).toHaveBeenCalled();
		const directRequests = fetchSpy.mock.calls.map(([url]) => {
			const outer = new URL(url);
			return new URL(
				outer.searchParams.get('url') || outer.href,
			);
		});
		expect(directRequests.length).toBeGreaterThan(20);
		expect(directRequests.some((url) => url.searchParams.has('mmseq'))).toBe(
			true,
		);
		expect(
			directRequests.some((url) => url.searchParams.has('concentration')),
		).toBe(true);

		for (const direct of directRequests) {
			expect(direct.searchParams.get('mg')).toBe('4.6');
			expect(direct.searchParams.get('mono')).toBe('17.2');
		}
	});
});
