import { jest } from '@jest/globals';
import { readFileSync } from 'node:fs';

import {
	buildTmRequestParams,
	calculateMeltingTempDifferences,
	calculateSnapbackTmWittwer,
	createStem,
	getOligoTm,
	getPrimerTm,
	useForwardPrimer,
} from '../dist/script.js';
import {
	renderDeltaTmTable,
	renderSnapbackPrimer,
} from '../src/js/pages/resultsRender.js';
import { renderStemDiagram } from '../src/js/pages/resultsStemDiagram.js';
import { readStoredTmConditions } from '../src/js/shared/tmConditions.js';
import { validateTmConditions } from '../src/js/shared/validators.js';

function directRequestUrl(requestUrl) {
	const outer = new URL(requestUrl, 'https://usnapback.test');
	const proxiedUrl = outer.searchParams.get('url');
	return proxiedUrl ? new URL(proxiedUrl) : outer;
}

function successfulTmResponse(tm = 42.5) {
	return {
		ok: true,
		text: async () => `<html><body><tm>${tm}</tm></body></html>`,
	};
}

function storageWith(values = {}) {
	return {
		getItem(key) {
			return Object.prototype.hasOwnProperty.call(values, key)
				? String(values[key])
				: null;
		},
	};
}

function expectIonParams(url, magnesiumMm, monovalentMm) {
	expect(url.searchParams.get('mg')).toBe(String(magnesiumMm));
	expect(url.searchParams.get('mono')).toBe(String(monovalentMm));
	expect(url.searchParams.has('dntp')).toBe(false);
	expect(url.searchParams.has('dntpconc')).toBe(false);
}

function expectPrimaryOligoParams(url, magnesiumMm, monovalentMm) {
	expectIonParams(url, magnesiumMm, monovalentMm);
	expect(url.searchParams.get('otype')).toBe('oligo');
	expect(url.searchParams.get('concentration')).toBe('0.5');
	expect(url.searchParams.get('limitingconc')).toBe('0.5');
	expect(url.searchParams.get('saltcalctype')).toBe('bpdenominator');
}

describe('2026 result annotations', () => {
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

describe('2026 page copy', () => {
	test('prefills the desired wild-type Tm with 65 on first visit', () => {
		const html = readFileSync(
			new URL('../src/pages/desiredTm.html', import.meta.url),
			'utf8',
		);
		const page = new DOMParser().parseFromString(html, 'text/html');
		const desiredTm = page.getElementById('desiredTm');

		expect(desiredTm.value).toBe('65');
	});

	test('labels the desired Tm as wild-type and defines free magnesium', () => {
		const html = readFileSync(
			new URL('../src/pages/desiredTm.html', import.meta.url),
			'utf8',
		);
		const page = new DOMParser().parseFromString(html, 'text/html');
		const desiredTmLabel = page.querySelector('label[for="desiredTm"]');
		const magnesiumLabel = page.querySelector('label[for="magnesiumMm"]');
		const instructions = page
			.querySelector('.instructions')
			.textContent.replace(/\s+/g, ' ')
			.trim();

		expect(desiredTmLabel.textContent.replace(/\s+/g, ' ').trim()).toBe(
			'Desired wild-type Tm (°C)',
		);

		expect(magnesiumLabel.textContent.replace(/\s+/g, ' ').trim()).toBe(
			'Free Mg2+ (mM)',
		);
		expect(instructions).toMatch(
			/free Mg2\+ concentration \(total Mg2\+ − total dNTPs\)/i,
		);
		expect(page.body.textContent).not.toMatch(/unbound/i);
		expect(page.querySelectorAll('.field-help')).toHaveLength(0);
		expect(
			page.getElementById('magnesiumMm').hasAttribute('aria-describedby'),
		).toBe(false);
		expect(
			page.getElementById('monovalentMm').hasAttribute('aria-describedby'),
		).toBe(false);
	});

	test('mentions both ion inputs in the fourth start-page step', () => {
		const html = readFileSync(
			new URL('../src/pages/start.html', import.meta.url),
			'utf8',
		);
		const page = new DOMParser().parseFromString(html, 'text/html');
		const step = page.querySelectorAll('.start-description ol > li')[3];
		const copy = step.textContent.replace(/\s+/g, ' ').trim();

		expect(copy).toMatch(/free Mg2\+/i);
		expect(copy).toMatch(/total monovalent cations/i);
		expect(copy.length).toBeLessThan(150);
	});

	test('uses the requested concise results labels and temporary note', () => {
		const html = readFileSync(
			new URL('../src/pages/results.html', import.meta.url),
			'utf8',
		);
		const page = new DOMParser().parseFromString(html, 'text/html');
		const copy = page.body.textContent.replace(/\s+/g, ' ').trim();
		const caption = page
			.getElementById('deltaCaption')
			.textContent.replace(/\s+/g, ' ')
			.trim();

		expect(copy).toMatch(/Wittwer\/empirical comparison:/);
		expect(copy).not.toMatch(/Carl\/Wittwer/i);
		expect(page.getElementById('freeMagnesiumMm')).toBeNull();
		expect(page.getElementById('totalMonovalentMm')).toBeNull();
		const bootstrapScript = [...page.querySelectorAll('script[type="module"]')]
			.find((script) => script.textContent.includes('bootstrap_started'));
		expect(bootstrapScript).toBeDefined();
		expect(bootstrapScript.hasAttribute('src')).toBe(false);
		expect(bootstrapScript.textContent).toMatch(
			/resultsDiagnostics\.js\?v=\$\{releaseId\}/,
		);
		expect(bootstrapScript.textContent).toMatch(
			/results\.js\?v=\$\{releaseId\}/,
		);
		expect(bootstrapScript.textContent).toMatch(/module_import_failed/);
		expect(bootstrapScript.textContent).toMatch(/initialization_failed/);
		expect(
			page.querySelector('link[href*="styles/pages/results.css?v="]'),
		).not.toBeNull();
		expect(caption).toBe(
			"Note this is now calculated entirely for each option including the logic, so stem lengths may vary. It won't be included in the final program.",
		);
	});
});

describe('Tm condition inputs and requests', () => {
	afterEach(() => {
		jest.restoreAllMocks();
		sessionStorage.clear();
		document.body.innerHTML = '';
	});
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

	test('builds API parameters from explicit cation values', () => {
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

	test('reads saved ion choices and falls back to the requested defaults', () => {
		expect(readStoredTmConditions(storageWith())).toEqual({
			magnesiumMm: 3,
			monovalentMm: 13.7,
		});
		expect(
			readStoredTmConditions(
				storageWith({ magnesiumMm: 4.6, monovalentMm: 17.2 }),
			),
		).toEqual({ magnesiumMm: 4.6, monovalentMm: 17.2 });
		expect(
			readStoredTmConditions(
				storageWith({ magnesiumMm: '', monovalentMm: 'invalid' }),
			),
		).toEqual({ magnesiumMm: 3, monovalentMm: 13.7 });
	});

	test('sends the complete primary oligo payload', async () => {
		const fetchSpy = jest
			.spyOn(global, 'fetch')
			.mockResolvedValue(successfulTmResponse());

		await expect(
			getOligoTm('ACGTACGT', null, {
				magnesiumMm: 5.1,
				monovalentMm: 9.4,
			}),
		).resolves.toBe(42.5);

		const direct = directRequestUrl(fetchSpy.mock.calls[0][0]);
		expectPrimaryOligoParams(direct, 5.1, 9.4);
	});

	test('keeps the complete oligo payload on mismatch requests', async () => {
		const fetchSpy = jest.spyOn(global, 'fetch').mockResolvedValue({
			ok: true,
			text: async () =>
				'<html><body><tm>42.5</tm><mmtm>31.2</mmtm></body></html>',
		});

		await expect(
			getOligoTm(
				'ACGTACGT',
				{ position: 3, type: 'C' },
				{ magnesiumMm: 5.1, monovalentMm: 9.4 },
			),
		).resolves.toBe(31.2);

		const direct = directRequestUrl(fetchSpy.mock.calls[0][0]);
		expectPrimaryOligoParams(direct, 5.1, 9.4);
		expect(direct.searchParams.get('mmseq')).toBe('ACGGACGT');
	});

	test('passes selected ions through the Wittwer calculation', async () => {
		const fetchSpy = jest
			.spyOn(global, 'fetch')
			.mockResolvedValue(successfulTmResponse());

		await calculateSnapbackTmWittwer('ACGTACGT', 6, undefined, {
			magnesiumMm: 4.8,
			monovalentMm: 12.6,
		});

		const direct = directRequestUrl(fetchSpy.mock.calls[0][0]);
		expectPrimaryOligoParams(direct, 4.8, 12.6);
	});

	test('sends primer mode and ions without oligo concentrations', async () => {
		const fetchSpy = jest
			.spyOn(global, 'fetch')
			.mockResolvedValue(successfulTmResponse(61.25));

		await expect(
			getPrimerTm('ACGTACGT', {
				magnesiumMm: 4.4,
				monovalentMm: 16.8,
			}),
		).resolves.toBe(61.25);

		const direct = directRequestUrl(fetchSpy.mock.calls[0][0]);
		expectIonParams(direct, 4.4, 16.8);
		expect(direct.searchParams.get('otype')).toBe('primer');
		expect(direct.searchParams.get('saltcalctype')).toBe('bpdenominator');
		expect(direct.searchParams.has('concentration')).toBe(false);
		expect(direct.searchParams.has('limitingconc')).toBe(false);
	});

	test('uses the default ions for primer requests before choices are saved', async () => {
		const fetchSpy = jest
			.spyOn(global, 'fetch')
			.mockResolvedValue(successfulTmResponse(61.25));

		await getPrimerTm('ACGTACGT');

		const direct = directRequestUrl(fetchSpy.mock.calls[0][0]);
		expectIonParams(direct, 3, 13.7);
		expect(direct.searchParams.get('otype')).toBe('primer');
		expect(direct.searchParams.has('concentration')).toBe(false);
		expect(direct.searchParams.has('limitingconc')).toBe(false);
	});

	test('passes saved ion choices from both active primer-page callers', async () => {
		const html = readFileSync(
			new URL('../src/pages/primers.html', import.meta.url),
			'utf8',
		);
		const page = new DOMParser().parseFromString(html, 'text/html');
		document.body.innerHTML = page.body.innerHTML;
		sessionStorage.clear();
		sessionStorage.setItem('ampliconSeq', 'ACGT'.repeat(20));
		sessionStorage.setItem('primerFwdStart', '0');
		sessionStorage.setItem('primerFwdEnd', '19');
		sessionStorage.setItem('primerRevStart', '60');
		sessionStorage.setItem('primerRevEnd', '79');
		sessionStorage.setItem('magnesiumMm', '4.9');
		sessionStorage.setItem('monovalentMm', '18.3');
		const fetchSpy = jest
			.spyOn(global, 'fetch')
			.mockResolvedValue(successfulTmResponse(61.25));

		await import('../dist/js/pages/primers.js');
		await new Promise((resolve) => setTimeout(resolve, 0));

		expect(fetchSpy).toHaveBeenCalledTimes(2);
		for (const [requestUrl] of fetchSpy.mock.calls) {
			const direct = directRequestUrl(requestUrl);
			expectIonParams(direct, 4.9, 18.3);
			expect(direct.searchParams.get('otype')).toBe('primer');
			expect(direct.searchParams.has('concentration')).toBe(false);
			expect(direct.searchParams.has('limitingconc')).toBe(false);
		}

		sessionStorage.clear();
		document.body.innerHTML = '';
	});

	test('sends selected conditions through every primary Wittwer caller', async () => {
		const fetchSpy = jest.spyOn(global, 'fetch').mockImplementation((url) => {
			const direct = directRequestUrl(url);
			const seq = direct.searchParams.get('seq') || '';
			const tm = 20 + seq.length * 2;
			const mismatchTm = tm - 10;
			return Promise.resolve({
				ok: true,
				text: async () =>
					`<html><body><tm>${tm}</tm><mmtm>${mismatchTm}</mmtm></body></html>`,
			});
		});

		const sequence = 'ACGT'.repeat(25);
		const snvSite = { index: 50, variantBase: 'A' };
		const conditions = { magnesiumMm: 4.6, monovalentMm: 17.2 };

		await useForwardPrimer(sequence, snvSite, conditions);
		const { bestStemLoc } = await createStem(
			sequence,
			snvSite,
			{ primerLen: 20, compPrimerLen: 20 },
			'C',
			true,
			60,
			conditions,
		);
		await calculateMeltingTempDifferences(
			sequence,
			snvSite,
			bestStemLoc,
			true,
			conditions,
		);

		expect(fetchSpy).toHaveBeenCalled();
		const directRequests = fetchSpy.mock.calls.map(([url]) =>
			directRequestUrl(url),
		);
		expect(directRequests.length).toBeGreaterThan(20);
		expect(directRequests.some((url) => url.searchParams.has('mmseq'))).toBe(
			true,
		);

		for (const direct of directRequests) {
			expectPrimaryOligoParams(direct, 4.6, 17.2);
		}
	});
});
