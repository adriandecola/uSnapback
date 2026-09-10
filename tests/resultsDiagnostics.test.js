import { jest } from '@jest/globals';
import { readFileSync } from 'node:fs';

import {
	RESULTS_DIAGNOSTICS_RELEASE,
	RESULTS_DIAGNOSTICS_PREFIX,
	RESULTS_DIAGNOSTICS_STORAGE_KEY,
	createResultsDiagnostics,
	describeDiagnosticError,
	detectBrowserSummary,
	sanitizeDiagnosticData,
	sanitizeDiagnosticText,
	summarizeResultsInputs,
} from '../src/js/pages/resultsDiagnostics.js';

describe('results diagnostics', () => {
	test('emits stable one-line JSON records with timing and correlation fields', () => {
		let now = 100;
		const consoleTarget = {
			info: jest.fn(),
			warn: jest.fn(),
			error: jest.fn(),
		};
		const storage = { setItem: jest.fn() };
		const diagnostics = createResultsDiagnostics({
			releaseId: 'test-release',
			runId: 'test-run',
			clock: () => now,
			timestamp: () => '2026-09-09T12:00:00.000Z',
			consoleTarget,
			storage,
		});

		now = 112.34;
		diagnostics.info('calculation_started', { workerAvailable: true });

		expect(consoleTarget.info).toHaveBeenCalledTimes(1);
		const line = consoleTarget.info.mock.calls[0][0];
		expect(line.startsWith(`${RESULTS_DIAGNOSTICS_PREFIX} `)).toBe(true);
		const record = JSON.parse(line.slice(RESULTS_DIAGNOSTICS_PREFIX.length + 1));
		expect(record).toMatchObject({
			schemaVersion: 1,
			releaseId: 'test-release',
			runId: 'test-run',
			sequenceNumber: 1,
			timestamp: '2026-09-09T12:00:00.000Z',
			elapsedMs: 12.3,
			page: 'results',
			level: 'info',
			event: 'calculation_started',
			workerAvailable: true,
		});
		expect(storage.setItem).toHaveBeenCalledWith(
			RESULTS_DIAGNOSTICS_STORAGE_KEY,
			expect.any(String),
		);
	});

	test('redacts DNA, raw payload fields, and URL query strings', () => {
		const sequence = 'ACGTACGTACGTACGT';
		const error = new Error(
			`Could not process ${sequence} at https://example.test/results?sample=${sequence}`,
		);
		const details = sanitizeDiagnosticData({
			sequence,
			sequenceLength: sequence.length,
			result: { snapbackSeq: sequence },
			error: describeDiagnosticError(error),
		});
		const serialized = JSON.stringify(details);

		expect(serialized).not.toContain(sequence);
		expect(details.sequence).toBe('[redacted]');
		expect(details.sequenceLength).toBe(16);
		expect(details.result).toBe('[redacted]');
		expect(details.error.message).toContain('[DNA:16 bases]');
		expect(details.error.message).not.toContain('?sample=');
		expect(sanitizeDiagnosticText(sequence)).toBe('[DNA:16 bases]');
	});

	test('keeps short alleles and indexes out of every emitted report surface', () => {
		const consoleTarget = {
			info: jest.fn(),
			warn: jest.fn(),
			error: jest.fn(),
		};
		const storage = { setItem: jest.fn() };
		const diagnostics = createResultsDiagnostics({
			runId: 'privacy-test',
			consoleTarget,
			storage,
		});
		const error = new Error(
			'Invalid DNA sequence: "ACGT"; SNV at index 3; received {"index":3,"variantBase":"T"} from https://example.test/results?sample=ACGT',
		);

		diagnostics.error('worker_calculation_failed', {
			error,
			result: { sequence: 'ACGT', variantBase: 'T' },
		});

		const surfaces = [
			consoleTarget.error.mock.calls[0][0],
			storage.setItem.mock.calls.at(-1)[1],
			diagnostics.getReport(),
		];
		for (const surface of surfaces) {
			expect(surface).not.toContain('"ACGT"');
			expect(surface).not.toContain('SNV at index 3');
			expect(surface).not.toContain('"index":3');
			expect(surface).not.toContain('"variantBase":"T"');
			expect(surface).not.toContain('?sample=');
			expect(surface).not.toContain('"result":{"sequence"');
		}
	});

	test('summarizes calculation inputs without including sequence or allele data', () => {
		const summary = summarizeResultsInputs(
			{
				seq: 'ACGTACGTACGT',
				fwdLen: 20,
				revLen: 21,
				snvIndex: 6,
				snvBase: 'T',
			},
			{
				tmC: 65,
				tmConditions: { magnesiumMm: 2.2, monovalentMm: 13.7 },
			},
		);

		expect(summary).toEqual({
			ampliconLength: 12,
			forwardPrimerLength: 20,
			reversePrimerLength: 21,
			hasSnvSelection: true,
			desiredTm: 65,
			magnesiumMm: 2.2,
			monovalentMm: 13.7,
		});
		expect(JSON.stringify(summary)).not.toMatch(/ACGT|snvIndex|snvBase/);
	});

	test('reports only a coarse browser name and major version', () => {
		expect(
			detectBrowserSummary(
				'Mozilla/5.0 Chrome/140.0.0.0 Safari/537.36 Edg/140.0.0.0',
			),
		).toBe('Edge 140');
		expect(detectBrowserSummary('custom browser')).toBe('Unknown');
	});

	test('identifies failed same-origin resources without logging URL queries', () => {
		const handlers = new Map();
		const windowTarget = {
			location: {
				href: 'https://example.test/pages/results.html',
				origin: 'https://example.test',
			},
			addEventListener: jest.fn((event, handler) => {
				handlers.set(event, handler);
			}),
			removeEventListener: jest.fn(),
		};
		const consoleTarget = {
			info: jest.fn(),
			warn: jest.fn(),
			error: jest.fn(),
		};
		const diagnostics = createResultsDiagnostics({
			consoleTarget,
			storage: null,
		});
		diagnostics.installGlobalHandlers(windowTarget);

		handlers.get('error')({
			target: {
				tagName: 'SCRIPT',
				src: 'https://example.test/js/pages/results.js?v=private-value',
			},
		});

		const line = consoleTarget.error.mock.calls[0][0];
		const record = JSON.parse(
			line.slice(RESULTS_DIAGNOSTICS_PREFIX.length + 1),
		);
		expect(record).toMatchObject({
			event: 'resource_load_failed',
			element: 'SCRIPT',
			resourcePath: '/js/pages/results.js',
		});
		expect(line).not.toContain('private-value');
	});

	test('keeps the diagnostics release aligned with every changed cached asset', () => {
		const resultsHtml = readFileSync(
			new URL('../src/pages/results.html', import.meta.url),
			'utf8',
		);
		const resultsJs = readFileSync(
			new URL('../src/js/pages/results.js', import.meta.url),
			'utf8',
		);
		const workerJs = readFileSync(
			new URL('../src/js/workers/snapbackWorker.js', import.meta.url),
			'utf8',
		);
		const calculationJs = readFileSync(
			new URL('../src/js/pages/resultsCalculation.js', import.meta.url),
			'utf8',
		);

		expect(resultsHtml).toContain(
			`const releaseId = '${RESULTS_DIAGNOSTICS_RELEASE}'`,
		);
		expect(resultsHtml).toContain(
			`resultsDiagnostics.js?v=\${releaseId}`,
		);
		expect(resultsHtml).toContain(`results.js?v=\${releaseId}`);
		expect(resultsHtml).toContain('void (async () => {');
		expect(resultsHtml).toContain("'module_import_failed'");
		expect(resultsHtml).toContain(
			`results.css?v=${RESULTS_DIAGNOSTICS_RELEASE}`,
		);
		expect(resultsJs).toContain(
			`script.js?v=${RESULTS_DIAGNOSTICS_RELEASE}`,
		);
		expect(resultsJs).toContain(
			`resultsCalculation.js?v=${RESULTS_DIAGNOSTICS_RELEASE}`,
		);
		expect(resultsJs).toContain(
			`resultsDiagnostics.js?v=${RESULTS_DIAGNOSTICS_RELEASE}`,
		);
		expect(calculationJs).toContain(
			`resultsDiagnostics.js?v=${RESULTS_DIAGNOSTICS_RELEASE}`,
		);
		expect(calculationJs).toContain(
			`snapbackWorkerProtocol.js?v=${RESULTS_DIAGNOSTICS_RELEASE}`,
		);
		expect(calculationJs).toContain(
			'`../workers/snapbackWorker.js?v=${RESULTS_DIAGNOSTICS_RELEASE}`',
		);
		expect(workerJs).toContain(
			`script.js?v=${RESULTS_DIAGNOSTICS_RELEASE}`,
		);
		expect(workerJs).toContain(
			`snapbackWorkerProtocol.js?v=${RESULTS_DIAGNOSTICS_RELEASE}`,
		);
	});

	test('does not retain legacy console dumps containing assay sequences', () => {
		const calculationSource = readFileSync(
			new URL('../src/script.js', import.meta.url),
			'utf8',
		);
		const resultsSource = readFileSync(
			new URL('../src/js/pages/results.js', import.meta.url),
			'utf8',
		);

		expect(calculationSource).not.toMatch(
			/console\.(?:log|debug|table|group|groupEnd)/,
		);
		expect(calculationSource).not.toMatch(
			/Stem \(WILD\) sequence|PARSING THERMO PARAMS RESPONSE HTML/,
		);
		expect(resultsSource).not.toMatch(/console\.log\(result\)/);
	});
});
