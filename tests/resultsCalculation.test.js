import { jest } from '@jest/globals';

import {
	calculateSnapbackForResults,
	waitForInitialPaint,
} from '../src/js/pages/resultsCalculation.js';
import {
	SNAPBACK_WORKER_CALCULATE,
	SNAPBACK_WORKER_CALCULATION_STARTED,
	SNAPBACK_WORKER_FAILURE,
	SNAPBACK_WORKER_READY,
	SNAPBACK_WORKER_SUCCESS,
} from '../src/js/shared/snapbackWorkerProtocol.js';

class FakeWorker {
	static latest = null;

	constructor(url, options) {
		this.url = url;
		this.options = options;
		this.onmessage = null;
		this.onerror = null;
		this.onmessageerror = null;
		this.terminated = false;
		this.sent = [];
		FakeWorker.latest = this;
	}

	postMessage(message) {
		this.sent.push(message);
	}

	terminate() {
		this.terminated = true;
	}
}

function createDiagnosticRecorder() {
	const records = [];
	const record = (level) => (event, details = {}) => {
		records.push({ level, event, ...details });
	};
	return {
		records,
		diagnostics: {
			info: record('info'),
			warn: record('warn'),
			error: record('error'),
			describeError: (error) => ({
				name: error?.name || 'Error',
				message: error?.message || String(error),
				...(error?.code ? { code: error.code } : {}),
			}),
		},
	};
}

describe('results calculation orchestration', () => {
	beforeEach(() => {
		jest.spyOn(console, 'info').mockImplementation(() => {});
		jest.spyOn(console, 'warn').mockImplementation(() => {});
		jest.spyOn(console, 'error').mockImplementation(() => {});
	});

	afterEach(() => {
		jest.useRealTimers();
		jest.restoreAllMocks();
		FakeWorker.latest = null;
	});

	test('finishes the initial-paint wait when requestAnimationFrame stays paused', async () => {
		jest.useFakeTimers();
		const requestFrame = jest.fn();
		const { diagnostics, records } = createDiagnosticRecorder();
		const waiting = waitForInitialPaint({
			requestFrame,
			timeoutMs: 25,
			diagnostics,
		});

		expect(requestFrame).toHaveBeenCalledTimes(1);
		jest.advanceTimersByTime(24);
		let finished = false;
		waiting.then(() => {
			finished = true;
		});
		await Promise.resolve();
		expect(finished).toBe(false);

		jest.advanceTimersByTime(1);
		await expect(waiting).resolves.toBeUndefined();
		expect(records).toEqual([
			expect.objectContaining({
				event: 'initial_paint_completed',
				via: 'timeout',
			}),
		]);
	});

	test('runs createSnapback in a module worker and terminates it on success', async () => {
		const directCalculate = jest.fn();
		const args = ['ACGT', 12, 12, { index: 2, variantBase: 'A' }, 60, {}];
		const { diagnostics, records } = createDiagnosticRecorder();
		let now = 10;
		const pending = calculateSnapbackForResults(args, {
			directCalculate,
			WorkerConstructor: FakeWorker,
			diagnostics,
			clock: () => now,
		});
		const worker = FakeWorker.latest;

		expect(worker.options).toEqual({
			type: 'module',
			name: 'usnapback-calculation',
		});
		expect(worker.sent).toEqual([
			{ type: SNAPBACK_WORKER_CALCULATE, args },
		]);

		now = 12;
		worker.onmessage({ data: { type: SNAPBACK_WORKER_READY } });
		now = 13;
		worker.onmessage({
			data: { type: SNAPBACK_WORKER_CALCULATION_STARTED },
		});
		now = 18;
		worker.onmessage({
			data: { type: SNAPBACK_WORKER_SUCCESS, result: { snapbackSeq: 'AC' } },
		});

		await expect(pending).resolves.toEqual({ snapbackSeq: 'AC' });
		expect(directCalculate).not.toHaveBeenCalled();
		expect(worker.terminated).toBe(true);
		expect(records.map(({ event }) => event)).toEqual([
			'calculation_started',
			'worker_created',
			'worker_request_sent',
			'worker_ready',
			'worker_calculation_started',
			'worker_succeeded',
			'calculation_succeeded',
		]);
		expect(records.at(-1)).toMatchObject({
			event: 'calculation_succeeded',
			mode: 'worker',
			durationMs: 8,
		});
		expect(records.filter(({ level }) => level === 'error')).toEqual([]);
		expect(JSON.stringify(records)).not.toContain('ACGT');
	});

	test('preserves a calculation error sent by the worker', async () => {
		const directCalculate = jest.fn();
		const pending = calculateSnapbackForResults([], {
			directCalculate,
			WorkerConstructor: FakeWorker,
		});
		const worker = FakeWorker.latest;

		worker.onmessage({
			data: {
				type: SNAPBACK_WORKER_FAILURE,
				error: {
					name: 'RangeError',
					message: 'No admissible stem.',
					code: 'NO_ADMISSIBLE_STEM',
					highestWildTm: 37.5,
				},
			},
		});

		await expect(pending).rejects.toMatchObject({
			name: 'RangeError',
			message: 'No admissible stem.',
			code: 'NO_ADMISSIBLE_STEM',
			highestWildTm: 37.5,
		});
		expect(directCalculate).not.toHaveBeenCalled();
		expect(worker.terminated).toBe(true);
	});

	test('uses the direct implementation when Worker is unavailable', async () => {
		const directCalculate = jest.fn().mockResolvedValue({ ok: true });
		const { diagnostics, records } = createDiagnosticRecorder();

		await expect(
			calculateSnapbackForResults(['sequence'], {
				directCalculate,
				WorkerConstructor: null,
				diagnostics,
			}),
		).resolves.toEqual({ ok: true });
		expect(directCalculate).toHaveBeenCalledWith('sequence');
		expect(records).toEqual(
			expect.arrayContaining([
				expect.objectContaining({
					event: 'direct_fallback_started',
					reason: 'worker_unavailable',
				}),
				expect.objectContaining({
					event: 'calculation_succeeded',
					mode: 'direct',
				}),
			]),
		);
		expect(records.at(-1)).toMatchObject({
			event: 'calculation_succeeded',
			mode: 'direct',
		});
		expect(records.filter(({ level }) => level === 'error')).toEqual([]);
	});

	test('falls back directly when a worker fails to load', async () => {
		const directCalculate = jest.fn().mockResolvedValue({ ok: true });
		const { diagnostics, records } = createDiagnosticRecorder();
		const pending = calculateSnapbackForResults(['sequence'], {
			directCalculate,
			WorkerConstructor: FakeWorker,
			diagnostics,
		});
		const worker = FakeWorker.latest;

		worker.onerror({
			message: 'Blocked by Content Security Policy',
			preventDefault: jest.fn(),
		});

		await expect(pending).resolves.toEqual({ ok: true });
		expect(directCalculate).toHaveBeenCalledWith('sequence');
		expect(worker.terminated).toBe(true);
		expect(records).toEqual(
			expect.arrayContaining([
				expect.objectContaining({
					event: 'worker_failed',
					reason: 'worker_load_or_runtime_error',
					error: expect.objectContaining({
						message: 'Blocked by Content Security Policy',
					}),
				}),
				expect.objectContaining({
					event: 'calculation_succeeded',
					mode: 'direct',
				}),
			]),
		);
	});

	test('terminates and rejects a worker that never responds', async () => {
		jest.useFakeTimers();
		const directCalculate = jest.fn();
		const { diagnostics, records } = createDiagnosticRecorder();
		const pending = calculateSnapbackForResults(['sequence'], {
			directCalculate,
			WorkerConstructor: FakeWorker,
			workerTimeoutMs: 60_000,
			diagnostics,
		});
		const worker = FakeWorker.latest;
		const rejection = expect(pending).rejects.toMatchObject({
			code: 'SNAPBACK_CALCULATION_TIMEOUT',
			message: expect.stringMatching(/retry.*shorten the amplicon/i),
		});

		jest.advanceTimersByTime(59_999);
		expect(worker.terminated).toBe(false);
		jest.advanceTimersByTime(1);

		await rejection;
		expect(worker.terminated).toBe(true);
		expect(directCalculate).not.toHaveBeenCalled();
		expect(jest.getTimerCount()).toBe(0);
		expect(records).toEqual(
			expect.arrayContaining([
				expect.objectContaining({
					event: 'worker_timeout',
					timeoutMs: 60_000,
					workerReady: false,
					calculationAcknowledged: false,
					error: expect.objectContaining({
						code: 'SNAPBACK_CALCULATION_TIMEOUT',
					}),
				}),
				expect.objectContaining({
					event: 'calculation_failed',
					mode: 'worker',
				}),
			]),
		);
		expect(
			records.some(({ event }) =>
				['calculation_succeeded', 'direct_fallback_started'].includes(event),
			),
		).toBe(false);
	});

	test('distinguishes a calculation stall from a worker startup stall', async () => {
		jest.useFakeTimers();
		const { diagnostics, records } = createDiagnosticRecorder();
		const pending = calculateSnapbackForResults(['sequence'], {
			directCalculate: jest.fn(),
			WorkerConstructor: FakeWorker,
			workerTimeoutMs: 25,
			diagnostics,
		});
		const worker = FakeWorker.latest;
		const rejection = expect(pending).rejects.toMatchObject({
			code: 'SNAPBACK_CALCULATION_TIMEOUT',
		});

		worker.onmessage({ data: { type: SNAPBACK_WORKER_READY } });
		worker.onmessage({
			data: { type: SNAPBACK_WORKER_CALCULATION_STARTED },
		});
		jest.advanceTimersByTime(25);

		await rejection;
		expect(records.find(({ event }) => event === 'worker_timeout')).toMatchObject(
			{
				workerReady: true,
				calculationAcknowledged: true,
			},
		);
	});
});
