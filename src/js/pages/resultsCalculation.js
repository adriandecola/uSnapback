import {
	SNAPBACK_WORKER_CALCULATE,
	SNAPBACK_WORKER_FAILURE,
	SNAPBACK_WORKER_SUCCESS,
	deserializeSnapbackWorkerError,
} from '../shared/snapbackWorkerProtocol.js';

const INITIAL_PAINT_FALLBACK_MS = 100;
const WORKER_TIMEOUT_MS = 60_000;

export function waitForInitialPaint({
	requestFrame =
		typeof globalThis.requestAnimationFrame === 'function'
			? globalThis.requestAnimationFrame.bind(globalThis)
			: null,
	scheduleTimeout = globalThis.setTimeout.bind(globalThis),
	cancelTimeout = globalThis.clearTimeout.bind(globalThis),
	timeoutMs = INITIAL_PAINT_FALLBACK_MS,
} = {}) {
	return new Promise((resolve) => {
		let settled = false;
		let timeoutId;
		const finish = () => {
			if (settled) return;
			settled = true;
			if (timeoutId !== undefined) cancelTimeout(timeoutId);
			resolve();
		};

		// requestAnimationFrame may remain paused while a tab or iframe is hidden.
		// Keep a timer in flight so loading can never depend on that callback alone.
		timeoutId = scheduleTimeout(finish, requestFrame ? timeoutMs : 0);
		if (requestFrame) {
			try {
				requestFrame(finish);
			} catch {
				finish();
			}
		}
	});
}

function defaultWorkerConstructor() {
	return typeof globalThis.Worker === 'function' ? globalThis.Worker : null;
}

export async function calculateSnapbackForResults(
	args,
	{
		directCalculate,
		WorkerConstructor = defaultWorkerConstructor(),
		workerUrl = new URL('../workers/snapbackWorker.js', import.meta.url),
		workerTimeoutMs = WORKER_TIMEOUT_MS,
		scheduleTimeout = globalThis.setTimeout.bind(globalThis),
		cancelTimeout = globalThis.clearTimeout.bind(globalThis),
	} = {},
) {
	if (!Array.isArray(args)) {
		throw new TypeError('Snapback calculation arguments must be an array.');
	}
	if (typeof directCalculate !== 'function') {
		throw new TypeError('A direct snapback calculation fallback is required.');
	}

	if (typeof WorkerConstructor !== 'function') {
		return directCalculate(...args);
	}

	let worker;
	try {
		worker = new WorkerConstructor(workerUrl, {
			type: 'module',
			name: 'usnapback-calculation',
		});
	} catch (error) {
		console.warn(
			'Unable to start the snapback worker; using the direct calculation fallback.',
			error,
		);
		return directCalculate(...args);
	}

	return new Promise((resolve, reject) => {
		let settled = false;
		let workerTimeoutId;

		const cleanup = () => {
			if (workerTimeoutId !== undefined) {
				cancelTimeout(workerTimeoutId);
				workerTimeoutId = undefined;
			}
			worker.onmessage = null;
			worker.onerror = null;
			worker.onmessageerror = null;
			worker.terminate();
		};

		const finish = (callback, value) => {
			if (settled) return;
			settled = true;
			cleanup();
			callback(value);
		};

		const useDirectFallback = (reason) => {
			if (settled) return;
			settled = true;
			cleanup();
			console.warn(
				'Snapback worker became unavailable; using the direct calculation fallback.',
				reason,
			);
			Promise.resolve()
				.then(() => directCalculate(...args))
				.then(resolve, reject);
		};

		worker.onmessage = (event) => {
			const message = event?.data;
			if (message?.type === SNAPBACK_WORKER_SUCCESS) {
				finish(resolve, message.result);
				return;
			}
			if (message?.type === SNAPBACK_WORKER_FAILURE) {
				finish(reject, deserializeSnapbackWorkerError(message.error));
				return;
			}

			useDirectFallback(
				new Error('The snapback worker returned an invalid response.'),
			);
		};

		worker.onerror = (event) => {
			event?.preventDefault?.();
			useDirectFallback(
				event?.error ||
					new Error(event?.message || 'The snapback worker failed to load.'),
			);
		};

		worker.onmessageerror = () => {
			useDirectFallback(
				new Error('The snapback worker returned unreadable data.'),
			);
		};

		workerTimeoutId = scheduleTimeout(() => {
			const error = new Error(
				'Snapback calculation took too long. Please retry, or shorten the amplicon if the problem continues.',
			);
			error.code = 'SNAPBACK_CALCULATION_TIMEOUT';
			finish(reject, error);
		}, workerTimeoutMs);

		try {
			worker.postMessage({
				type: SNAPBACK_WORKER_CALCULATE,
				args,
			});
		} catch (error) {
			useDirectFallback(error);
		}
	});
}
