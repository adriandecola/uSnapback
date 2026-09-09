import { createSnapback } from '../../script.js';
import {
	SNAPBACK_WORKER_CALCULATE,
	SNAPBACK_WORKER_FAILURE,
	SNAPBACK_WORKER_SUCCESS,
	serializeSnapbackWorkerError,
} from '../shared/snapbackWorkerProtocol.js';

globalThis.addEventListener('message', async (event) => {
	const message = event?.data;
	if (message?.type !== SNAPBACK_WORKER_CALCULATE) return;

	try {
		if (!Array.isArray(message.args)) {
			throw new TypeError('Snapback calculation arguments must be an array.');
		}
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
