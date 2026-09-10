import {
	SNAPBACK_WORKER_CALCULATE,
	SNAPBACK_WORKER_CALCULATION_STARTED,
	SNAPBACK_WORKER_FAILURE,
	SNAPBACK_WORKER_READY,
	SNAPBACK_WORKER_SUCCESS,
	deserializeSnapbackWorkerError,
} from '../shared/snapbackWorkerProtocol.js?v=20260909.2';
import {
	RESULTS_DIAGNOSTICS_RELEASE,
	RESULTS_DIAGNOSTICS_PREFIX,
	describeDiagnosticError,
	sanitizeDiagnosticData,
} from './resultsDiagnostics.js?v=20260909.2';

const INITIAL_PAINT_FALLBACK_MS = 100;
const WORKER_TIMEOUT_MS = 60_000;

function defaultClock() {
	return typeof globalThis.performance?.now === 'function'
		? globalThis.performance.now()
		: Date.now();
}

function durationSince(clock, startedAt) {
	return Math.round((clock() - startedAt) * 10) / 10;
}

function reportDiagnostic(diagnostics, level, event, details = {}) {
	try {
		if (typeof diagnostics?.[level] === 'function') {
			diagnostics[level](event, details);
			return;
		}
		const writer = globalThis.console?.[level] || globalThis.console?.log;
		writer?.call(
			globalThis.console,
			`${RESULTS_DIAGNOSTICS_PREFIX} ${JSON.stringify({
				level,
				event,
				...sanitizeDiagnosticData(details),
			})}`,
		);
	} catch {
		// Diagnostics must never alter the calculation path.
	}
}

function diagnosticError(diagnostics, error) {
	return typeof diagnostics?.describeError === 'function'
		? diagnostics.describeError(error)
		: describeDiagnosticError(error);
}

export function waitForInitialPaint({
	requestFrame =
		typeof globalThis.requestAnimationFrame === 'function'
			? globalThis.requestAnimationFrame.bind(globalThis)
			: null,
	scheduleTimeout = globalThis.setTimeout.bind(globalThis),
	cancelTimeout = globalThis.clearTimeout.bind(globalThis),
	timeoutMs = INITIAL_PAINT_FALLBACK_MS,
	diagnostics = null,
	clock = defaultClock,
} = {}) {
	const startedAt = clock();
	return new Promise((resolve) => {
		let settled = false;
		let timeoutId;
		const finish = (via, error = null) => {
			if (settled) return;
			settled = true;
			if (timeoutId !== undefined) cancelTimeout(timeoutId);
			reportDiagnostic(diagnostics, 'info', 'initial_paint_completed', {
				via,
				durationMs: durationSince(clock, startedAt),
				...(error
					? { error: diagnosticError(diagnostics, error) }
					: {}),
			});
			resolve();
		};

		// requestAnimationFrame may remain paused while a tab or iframe is hidden.
		// Keep a timer in flight so loading can never depend on that callback alone.
		timeoutId = scheduleTimeout(
			() => finish(requestFrame ? 'timeout' : 'no_animation_frame'),
			requestFrame ? timeoutMs : 0,
		);
		if (requestFrame) {
			try {
				requestFrame(() => finish('animation_frame'));
			} catch (error) {
				finish('animation_frame_error', error);
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
		workerUrl = new URL(
			`../workers/snapbackWorker.js?v=${RESULTS_DIAGNOSTICS_RELEASE}`,
			import.meta.url,
		),
		workerTimeoutMs = WORKER_TIMEOUT_MS,
		scheduleTimeout = globalThis.setTimeout.bind(globalThis),
		cancelTimeout = globalThis.clearTimeout.bind(globalThis),
		diagnostics = null,
		clock = defaultClock,
	} = {},
) {
	if (!Array.isArray(args)) {
		throw new TypeError('Snapback calculation arguments must be an array.');
	}
	if (typeof directCalculate !== 'function') {
		throw new TypeError('A direct snapback calculation fallback is required.');
	}
	const calculationStartedAt = clock();
	reportDiagnostic(diagnostics, 'info', 'calculation_started', {
		workerAvailable: typeof WorkerConstructor === 'function',
		timeoutMs: workerTimeoutMs,
	});

	const runDirectCalculation = async (reason, workerError = null) => {
		const directStartedAt = clock();
		reportDiagnostic(diagnostics, 'warn', 'direct_fallback_started', {
			reason,
			...(workerError
				? { error: diagnosticError(diagnostics, workerError) }
				: {}),
		});
		try {
			const result = await directCalculate(...args);
			reportDiagnostic(diagnostics, 'info', 'direct_fallback_succeeded', {
				reason,
				durationMs: durationSince(clock, directStartedAt),
			});
			reportDiagnostic(diagnostics, 'info', 'calculation_succeeded', {
				mode: 'direct',
				durationMs: durationSince(clock, calculationStartedAt),
			});
			return result;
		} catch (error) {
			reportDiagnostic(diagnostics, 'error', 'direct_fallback_failed', {
				reason,
				durationMs: durationSince(clock, directStartedAt),
				error: diagnosticError(diagnostics, error),
			});
			reportDiagnostic(diagnostics, 'error', 'calculation_failed', {
				mode: 'direct',
				durationMs: durationSince(clock, calculationStartedAt),
				error: diagnosticError(diagnostics, error),
			});
			throw error;
		}
	};

	if (typeof WorkerConstructor !== 'function') {
		return runDirectCalculation('worker_unavailable');
	}

	let worker;
	try {
		worker = new WorkerConstructor(workerUrl, {
			type: 'module',
			name: 'usnapback-calculation',
		});
		reportDiagnostic(diagnostics, 'info', 'worker_created', {
			durationMs: durationSince(clock, calculationStartedAt),
		});
	} catch (error) {
		reportDiagnostic(diagnostics, 'warn', 'worker_creation_failed', {
			error: diagnosticError(diagnostics, error),
		});
		return runDirectCalculation('worker_creation_failed', error);
	}

	return new Promise((resolve, reject) => {
		let settled = false;
		let workerTimeoutId;
		let workerReady = false;
		let calculationAcknowledged = false;

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

		const useDirectFallback = (reason, error) => {
			if (settled) return;
			settled = true;
			cleanup();
			reportDiagnostic(diagnostics, 'warn', 'worker_failed', {
				reason,
				durationMs: durationSince(clock, calculationStartedAt),
				error: diagnosticError(diagnostics, error),
			});
			Promise.resolve()
				.then(() => runDirectCalculation(reason, error))
				.then(resolve, reject);
		};

		worker.onmessage = (event) => {
			const message = event?.data;
			if (message?.type === SNAPBACK_WORKER_READY) {
				workerReady = true;
				reportDiagnostic(diagnostics, 'info', 'worker_ready', {
					durationMs: durationSince(clock, calculationStartedAt),
				});
				return;
			}
			if (message?.type === SNAPBACK_WORKER_CALCULATION_STARTED) {
				workerReady = true;
				calculationAcknowledged = true;
				reportDiagnostic(diagnostics, 'info', 'worker_calculation_started', {
					durationMs: durationSince(clock, calculationStartedAt),
				});
				return;
			}
			if (message?.type === SNAPBACK_WORKER_SUCCESS) {
				reportDiagnostic(diagnostics, 'info', 'worker_succeeded', {
					durationMs: durationSince(clock, calculationStartedAt),
				});
				reportDiagnostic(diagnostics, 'info', 'calculation_succeeded', {
					mode: 'worker',
					durationMs: durationSince(clock, calculationStartedAt),
				});
				finish(resolve, message.result);
				return;
			}
			if (message?.type === SNAPBACK_WORKER_FAILURE) {
				const error = deserializeSnapbackWorkerError(message.error);
				reportDiagnostic(diagnostics, 'error', 'worker_calculation_failed', {
					durationMs: durationSince(clock, calculationStartedAt),
					error: diagnosticError(diagnostics, error),
				});
				reportDiagnostic(diagnostics, 'error', 'calculation_failed', {
					mode: 'worker',
					durationMs: durationSince(clock, calculationStartedAt),
					error: diagnosticError(diagnostics, error),
				});
				finish(reject, error);
				return;
			}

			useDirectFallback(
				'invalid_worker_response',
				new Error('The snapback worker returned an invalid response.'),
			);
		};

		worker.onerror = (event) => {
			event?.preventDefault?.();
			useDirectFallback(
				'worker_load_or_runtime_error',
				event?.error ||
					new Error(event?.message || 'The snapback worker failed to load.'),
			);
		};

		worker.onmessageerror = () => {
			useDirectFallback(
				'worker_message_error',
				new Error('The snapback worker returned unreadable data.'),
			);
		};

		workerTimeoutId = scheduleTimeout(() => {
			const error = new Error(
				'Snapback calculation took too long. Please retry, or shorten the amplicon if the problem continues.',
			);
			error.code = 'SNAPBACK_CALCULATION_TIMEOUT';
			reportDiagnostic(diagnostics, 'error', 'worker_timeout', {
				timeoutMs: workerTimeoutMs,
				durationMs: durationSince(clock, calculationStartedAt),
				workerReady,
				calculationAcknowledged,
				error: diagnosticError(diagnostics, error),
			});
			reportDiagnostic(diagnostics, 'error', 'calculation_failed', {
				mode: 'worker',
				durationMs: durationSince(clock, calculationStartedAt),
				error: diagnosticError(diagnostics, error),
			});
			finish(reject, error);
		}, workerTimeoutMs);

		try {
			worker.postMessage({
				type: SNAPBACK_WORKER_CALCULATE,
				args,
			});
			reportDiagnostic(diagnostics, 'info', 'worker_request_sent');
		} catch (error) {
			useDirectFallback('worker_request_failed', error);
		}
	});
}
