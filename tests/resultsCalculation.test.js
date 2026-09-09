import { jest } from '@jest/globals';

import {
	calculateSnapbackForResults,
	waitForInitialPaint,
} from '../src/js/pages/resultsCalculation.js';
import {
	SNAPBACK_WORKER_CALCULATE,
	SNAPBACK_WORKER_FAILURE,
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

describe('results calculation orchestration', () => {
	afterEach(() => {
		jest.useRealTimers();
		jest.restoreAllMocks();
		FakeWorker.latest = null;
	});

	test('finishes the initial-paint wait when requestAnimationFrame stays paused', async () => {
		jest.useFakeTimers();
		const requestFrame = jest.fn();
		const waiting = waitForInitialPaint({ requestFrame, timeoutMs: 25 });

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
	});

	test('runs createSnapback in a module worker and terminates it on success', async () => {
		const directCalculate = jest.fn();
		const args = ['ACGT', 12, 12, { index: 2, variantBase: 'A' }, 60, {}];
		const pending = calculateSnapbackForResults(args, {
			directCalculate,
			WorkerConstructor: FakeWorker,
		});
		const worker = FakeWorker.latest;

		expect(worker.options).toEqual({
			type: 'module',
			name: 'usnapback-calculation',
		});
		expect(worker.sent).toEqual([
			{ type: SNAPBACK_WORKER_CALCULATE, args },
		]);

		worker.onmessage({
			data: { type: SNAPBACK_WORKER_SUCCESS, result: { snapbackSeq: 'AC' } },
		});

		await expect(pending).resolves.toEqual({ snapbackSeq: 'AC' });
		expect(directCalculate).not.toHaveBeenCalled();
		expect(worker.terminated).toBe(true);
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

		await expect(
			calculateSnapbackForResults(['sequence'], {
				directCalculate,
				WorkerConstructor: null,
			}),
		).resolves.toEqual({ ok: true });
		expect(directCalculate).toHaveBeenCalledWith('sequence');
	});

	test('falls back directly when a worker fails to load', async () => {
		const warn = jest.spyOn(console, 'warn').mockImplementation(() => {});
		const directCalculate = jest.fn().mockResolvedValue({ ok: true });
		const pending = calculateSnapbackForResults(['sequence'], {
			directCalculate,
			WorkerConstructor: FakeWorker,
		});
		const worker = FakeWorker.latest;

		worker.onerror({
			message: 'Blocked by Content Security Policy',
			preventDefault: jest.fn(),
		});

		await expect(pending).resolves.toEqual({ ok: true });
		expect(directCalculate).toHaveBeenCalledWith('sequence');
		expect(worker.terminated).toBe(true);
		expect(warn).toHaveBeenCalledTimes(1);
	});

	test('terminates and rejects a worker that never responds', async () => {
		jest.useFakeTimers();
		const directCalculate = jest.fn();
		const pending = calculateSnapbackForResults(['sequence'], {
			directCalculate,
			WorkerConstructor: FakeWorker,
			workerTimeoutMs: 60_000,
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
	});
});
