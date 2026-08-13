import { jest } from '@jest/globals';

import {
	buildTmRequestParams,
	getOligoTm,
} from '../dist/script.js';
import {
	renderDeltaTmTable,
	renderSnapbackPrimer,
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
			document.querySelector('.seq-seg--inner-loop-mismatch').textContent,
		).toBe('C');
		expect(document.querySelector('.seq-seg--tail').textContent).toBe(
			'ACGT',
		);
		expect(document.querySelector('.seq-seg--primer').textContent).toBe(
			'TTAA',
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

	test('adds allele and tail bases to both axes of the delta-Tm table', () => {
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
			'Wild-type match (A)',
		);
		expect(document.getElementById('dt-var-heading').textContent).toBe(
			'Variant match (G)',
		);
		expect(document.getElementById('dt-fwd-heading').textContent).toBe(
			'Tail on forward primer (T/C)',
		);
		expect(document.getElementById('dt-rev-heading').textContent).toBe(
			'Tail on reverse primer (A/G)',
		);
	});
});

describe('Carl-requested cation inputs', () => {
	test('accepts the requested default free-magnesium and monovalent values', () => {
		expect(validateTmConditions('3.0', '13.7')).toEqual({
			ok: true,
			msg: '',
			data: { magnesiumMm: 3, monovalentMm: 13.7 },
		});
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
});
