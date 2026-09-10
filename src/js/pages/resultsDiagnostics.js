export const RESULTS_DIAGNOSTICS_RELEASE = '20260909.2';
export const RESULTS_DIAGNOSTICS_PREFIX = '[uSnapback diagnostics]';
export const RESULTS_DIAGNOSTICS_STORAGE_KEY =
	'usnapback.resultsDiagnostics';

const DIAGNOSTICS_SCHEMA_VERSION = 1;
const MAX_RECORDS = 100;
const MAX_TEXT_LENGTH = 4000;
const DNA_RUN_PATTERN = /[ACGT]{8,}/gi;
const QUOTED_DNA_PATTERN = /(["'`])([ACGT]+)\1/gi;
const SENSITIVE_JSON_NUMBER_PATTERN =
	/("(?:index|position)"\s*:\s*)-?\d+/gi;
const SENSITIVE_INDEX_PATTERN =
	/\b(SNV(?:\s+at)?\s+index|snvSite\.index|mismatchPos|mismatch(?:\.position|\s+position))(\s*(?:\(|:)?\s*)-?\d+/gi;
const URL_QUERY_PATTERN = /(https?:\/\/[^\s?#]+)[?#][^\s)]+/gi;

function defaultClock() {
	return typeof globalThis.performance?.now === 'function'
		? globalThis.performance.now()
		: Date.now();
}

function defaultTimestamp() {
	return new Date().toISOString();
}

function defaultStorage() {
	try {
		return globalThis.sessionStorage || null;
	} catch {
		return null;
	}
}

function createRunId() {
	try {
		if (typeof globalThis.crypto?.randomUUID === 'function') {
			return globalThis.crypto.randomUUID();
		}
	} catch {
		// Fall through to the non-cryptographic correlation ID below.
	}
	return `${Date.now().toString(36)}-${Math.random().toString(36).slice(2, 10)}`;
}

function roundMilliseconds(value) {
	return Math.round(value * 10) / 10;
}

export function sanitizeDiagnosticText(value) {
	const sanitized = String(value)
		.replace(
			QUOTED_DNA_PATTERN,
			(_match, _quote, sequence) => `[DNA:${sequence.length} bases]`,
		)
		.replace(DNA_RUN_PATTERN, (sequence) => `[DNA:${sequence.length} bases]`)
		.replace(SENSITIVE_JSON_NUMBER_PATTERN, '$1"[redacted]"')
		.replace(
			SENSITIVE_INDEX_PATTERN,
			(_match, label, separator) => `${label}${separator}[redacted]`,
		)
		.replace(URL_QUERY_PATTERN, '$1');
	return sanitized.length <= MAX_TEXT_LENGTH
		? sanitized
		: `${sanitized.slice(0, MAX_TEXT_LENGTH)}…`;
}

export function describeDiagnosticError(error) {
	const errorLike =
		error != null && (typeof error === 'object' || typeof error === 'function')
			? error
			: null;
	const details = {
		name: sanitizeDiagnosticText(errorLike?.name || 'Error'),
		message: sanitizeDiagnosticText(
			errorLike?.message ||
				(error == null ? 'Unknown error' : String(error)),
		),
	};
	if (errorLike?.code != null) {
		details.code = sanitizeDiagnosticText(errorLike.code);
	}
	if (errorLike?.stack) {
		details.stack = sanitizeDiagnosticText(errorLike.stack);
	}
	return details;
}

function shouldRedactKey(key) {
	const normalized = String(key).toLowerCase();
	if (normalized.endsWith('length')) return false;
	return (
		normalized === 'seq' ||
		normalized.includes('sequence') ||
		normalized.endsWith('base') ||
		normalized === 'result' ||
		normalized === 'payload' ||
		normalized === 'args' ||
		normalized === 'html' ||
		normalized === 'rawhtml'
	);
}

export function sanitizeDiagnosticData(value, depth = 0, seen = new WeakSet()) {
	if (value == null || typeof value === 'number' || typeof value === 'boolean') {
		return value;
	}
	if (typeof value === 'string') return sanitizeDiagnosticText(value);
	if (value instanceof Error) return describeDiagnosticError(value);
	if (typeof value !== 'object') return sanitizeDiagnosticText(value);
	if (seen.has(value)) return '[circular]';
	if (depth >= 4) return '[truncated]';

	seen.add(value);
	if (Array.isArray(value)) {
		return value
			.slice(0, 20)
			.map((item) => sanitizeDiagnosticData(item, depth + 1, seen));
	}

	const sanitized = {};
	for (const [key, item] of Object.entries(value).slice(0, 40)) {
		sanitized[key] = shouldRedactKey(key)
			? '[redacted]'
			: sanitizeDiagnosticData(item, depth + 1, seen);
	}
	return sanitized;
}

export function summarizeResultsInputs(inputs, validated) {
	return {
		ampliconLength:
			typeof inputs?.seq === 'string' ? inputs.seq.length : null,
		forwardPrimerLength: Number.isFinite(inputs?.fwdLen)
			? inputs.fwdLen
			: null,
		reversePrimerLength: Number.isFinite(inputs?.revLen)
			? inputs.revLen
			: null,
		hasSnvSelection:
			Number.isInteger(inputs?.snvIndex) &&
			typeof inputs?.snvBase === 'string' &&
			inputs.snvBase.length === 1,
		desiredTm: Number.isFinite(validated?.tmC) ? validated.tmC : null,
		magnesiumMm: Number.isFinite(validated?.tmConditions?.magnesiumMm)
			? validated.tmConditions.magnesiumMm
			: null,
		monovalentMm: Number.isFinite(validated?.tmConditions?.monovalentMm)
			? validated.tmConditions.monovalentMm
			: null,
	};
}

export function detectBrowserSummary(userAgent = '') {
	const ua = String(userAgent);
	for (const [name, pattern] of [
		['Edge', /Edg\/(\d+)/],
		['Chrome', /Chrome\/(\d+)/],
		['Firefox', /Firefox\/(\d+)/],
		['Safari', /Version\/(\d+).*Safari/],
	]) {
		const match = ua.match(pattern);
		if (match) return `${name} ${match[1]}`;
	}
	return 'Unknown';
}

function describeResourcePath(value, windowTarget) {
	if (!value) return null;
	try {
		const resourceUrl = new URL(value, windowTarget?.location?.href);
		const pageOrigin = windowTarget?.location?.origin;
		if (pageOrigin && resourceUrl.origin !== pageOrigin) {
			return '[cross-origin resource]';
		}
		return sanitizeDiagnosticText(resourceUrl.pathname || '/');
	} catch {
		return '[unavailable resource]';
	}
}

export function createResultsDiagnostics({
	releaseId = RESULTS_DIAGNOSTICS_RELEASE,
	runId = createRunId(),
	clock = defaultClock,
	timestamp = defaultTimestamp,
	consoleTarget = globalThis.console,
	storage = defaultStorage(),
} = {}) {
	const startedAt = clock();
	const records = [];
	let sequenceNumber = 0;
	let removeGlobalHandlers = null;

	const elapsedMs = () => roundMilliseconds(clock() - startedAt);

	const persist = () => {
		try {
			storage?.setItem(
				RESULTS_DIAGNOSTICS_STORAGE_KEY,
				JSON.stringify(records),
			);
		} catch {
			// Diagnostics must never interrupt the calculation if storage is blocked.
		}
	};

	const emit = (level, event, details = {}) => {
		sequenceNumber += 1;
		const sanitizedDetails = sanitizeDiagnosticData(details);
		for (const reservedKey of [
			'schemaVersion',
			'releaseId',
			'runId',
			'sequenceNumber',
			'timestamp',
			'elapsedMs',
			'page',
			'level',
			'event',
		]) {
			delete sanitizedDetails[reservedKey];
		}
		const record = {
			schemaVersion: DIAGNOSTICS_SCHEMA_VERSION,
			releaseId,
			runId,
			sequenceNumber,
			timestamp: timestamp(),
			elapsedMs: elapsedMs(),
			page: 'results',
			level,
			event: sanitizeDiagnosticText(event),
			...sanitizedDetails,
		};
		records.push(record);
		if (records.length > MAX_RECORDS) records.shift();
		persist();

		const writer = consoleTarget?.[level] || consoleTarget?.log;
		writer?.call(
			consoleTarget,
			`${RESULTS_DIAGNOSTICS_PREFIX} ${JSON.stringify(record)}`,
		);
		return record;
	};

	const diagnostics = {
		releaseId,
		runId,
		elapsedMs,
		describeError: describeDiagnosticError,
		info: (event, details) => emit('info', event, details),
		warn: (event, details) => emit('warn', event, details),
		error: (event, details) => emit('error', event, details),
		getReport: () => JSON.stringify(records, null, 2),
		installGlobalHandlers(windowTarget = globalThis.window) {
			if (!windowTarget?.addEventListener || removeGlobalHandlers) {
				return removeGlobalHandlers || (() => {});
			}

			const handleError = (event) => {
				if (event?.target && event.target !== windowTarget) {
					const resourcePath = describeResourcePath(
						event.target.currentSrc || event.target.src || event.target.href,
						windowTarget,
					);
					diagnostics.error('resource_load_failed', {
						element: event.target.tagName || 'unknown',
						...(resourcePath ? { resourcePath } : {}),
					});
					return;
				}
				const sourcePath = describeResourcePath(
					event?.filename,
					windowTarget,
				);
				diagnostics.error('window_error', {
					error: describeDiagnosticError(
						event?.error || event?.message || 'Unknown window error',
					),
					...(sourcePath ? { sourcePath } : {}),
					...(Number.isFinite(event?.lineno)
						? { lineNumber: event.lineno }
						: {}),
					...(Number.isFinite(event?.colno)
						? { columnNumber: event.colno }
						: {}),
				});
			};
			const handleRejection = (event) => {
				diagnostics.error('unhandled_rejection', {
					error: describeDiagnosticError(event?.reason),
				});
			};
			windowTarget.addEventListener('error', handleError, true);
			windowTarget.addEventListener(
				'unhandledrejection',
				handleRejection,
			);
			removeGlobalHandlers = () => {
				windowTarget.removeEventListener('error', handleError, true);
				windowTarget.removeEventListener(
					'unhandledrejection',
					handleRejection,
				);
				removeGlobalHandlers = null;
			};
			return removeGlobalHandlers;
		},
	};

	return diagnostics;
}
