export const SNAPBACK_WORKER_CALCULATE = 'usnapback:calculate';
export const SNAPBACK_WORKER_READY = 'usnapback:ready';
export const SNAPBACK_WORKER_CALCULATION_STARTED =
	'usnapback:calculation-started';
export const SNAPBACK_WORKER_SUCCESS = 'usnapback:success';
export const SNAPBACK_WORKER_FAILURE = 'usnapback:failure';

export function serializeSnapbackWorkerError(error) {
	const message =
		typeof error?.message === 'string' && error.message
			? error.message
			: String(error || 'Snapback calculation failed.');
	const payload = {
		name:
			typeof error?.name === 'string' && error.name ? error.name : 'Error',
		message,
	};

	if (typeof error?.stack === 'string') payload.stack = error.stack;
	if (typeof error?.code === 'string') payload.code = error.code;
	if (Number.isFinite(error?.highestWildTm)) {
		payload.highestWildTm = error.highestWildTm;
	}
	if (Array.isArray(error?.candidateFailures)) {
		payload.candidateFailures = error.candidateFailures;
	}

	return payload;
}

export function deserializeSnapbackWorkerError(payload) {
	const message =
		typeof payload?.message === 'string' && payload.message
			? payload.message
			: 'Snapback calculation failed in the background worker.';
	const error = new Error(message);

	if (typeof payload?.name === 'string' && payload.name) {
		error.name = payload.name;
	}
	if (typeof payload?.stack === 'string' && payload.stack) {
		error.stack = payload.stack;
	}
	if (typeof payload?.code === 'string') error.code = payload.code;
	if (Number.isFinite(payload?.highestWildTm)) {
		error.highestWildTm = payload.highestWildTm;
	}
	if (Array.isArray(payload?.candidateFailures)) {
		error.candidateFailures = payload.candidateFailures;
	}

	return error;
}
