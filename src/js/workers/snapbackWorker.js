import { createSnapback } from '../../script.js?v=20260909.2';
import {
	SNAPBACK_WORKER_CALCULATE,
	SNAPBACK_WORKER_CALCULATION_STARTED,
	SNAPBACK_WORKER_FAILURE,
	SNAPBACK_WORKER_READY,
	SNAPBACK_WORKER_SUCCESS,
	serializeSnapbackWorkerError,
} from '../shared/snapbackWorkerProtocol.js?v=20260909.2';

globalThis.postMessage({ type: SNAPBACK_WORKER_READY });

globalThis.addEventListener('message', async (event) => {
	const message = event?.data;
	if (message?.type !== SNAPBACK_WORKER_CALCULATE) return;

	try {
		if (!Array.isArray(message.args)) {
			throw new TypeError('Snapback calculation arguments must be an array.');
		}
		globalThis.postMessage({ type: SNAPBACK_WORKER_CALCULATION_STARTED });
		const result = await createSnapback(...message.args);
		globalThis.postMessage({
			type: SNAPBACK_WORKER_SUCCESS,
			result,
		});
	} catch (error) {
		globalThis.postMessage({
			type: SNAPBACK_WORKER_FAILURE,
			error: serializeSnapbackWorkerError(error),
		});
	}
});
