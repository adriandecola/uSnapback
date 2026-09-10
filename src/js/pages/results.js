/*
  File:             results.js
  Description:      Page logic for results.html
  Author:           Adrian deCola
  Relative Path:    uSnapback/src/js/pages/results.js
  Used by:          ../../pages/results.html
*/

/* ---------------------------------------- Imports --------------------------------------- */
import { createSnapback } from '../../script.js?v=20260909.2';
import {
	DEFAULT_MAGNESIUM_MM,
	DEFAULT_MONOVALENT_MM,
} from '../shared/constants.js';

import { wireCopyButton } from '../shared/clipboard.js';
import { renderStemDiagram } from './resultsStemDiagram.js';
import {
	calculateSnapbackForResults,
	waitForInitialPaint,
} from './resultsCalculation.js?v=20260909.2';
import {
	createResultsDiagnostics,
	describeDiagnosticError,
	detectBrowserSummary,
	summarizeResultsInputs,
} from './resultsDiagnostics.js?v=20260909.2';
import {
	renderSnapbackPrimer,
	renderLimitingPrimer,
	renderTailSummary,
	renderTmSummary,
	renderStemAndLoopSizes,
	renderDeltaTmTable,
} from './resultsRender.js';

// Importing validator functions
import {
	validateAmpliconSeq,
	validatePrimerLengths,
	validateSnv,
	validateDesiredTm,
	validateTmConditions,
} from '../shared/validators.js';

const PREV_PAGE = 'desiredTm.html';
const START_PAGE = 'start.html';
const AMPLICON_PAGE = 'amplicon.html';

/* --------------------------------------------------
Helper: recover to amplicon page when upstream inputs are invalid
-------------------------------------------------- */
function goBack(msg) {
	alert(
		msg +
			'\n\nYou will now be redirected to the amplicon page so you can adjust your inputs.',
	);
	window.location.href = AMPLICON_PAGE;
}

/* --------------------------------------------------
Helper: read inputs passed forward via sessionStorage
-------------------------------------------------- */
function readInputs() {
	return {
		seq: sessionStorage.getItem('ampliconSeqCropped') || '',
		fwdLen: +sessionStorage.getItem('forwardPrimerLen'),
		revLen: +sessionStorage.getItem('reversePrimerLen'),
		snvIndex: +sessionStorage.getItem('snvIndex'),
		snvBase: sessionStorage.getItem('snvBase'),
		tmStr: sessionStorage.getItem('desiredTm') ?? '',
		magnesiumStr:
			sessionStorage.getItem('magnesiumMm') ??
			String(DEFAULT_MAGNESIUM_MM),
		monovalentStr:
			sessionStorage.getItem('monovalentMm') ??
			String(DEFAULT_MONOVALENT_MM),
	};
}

/* --------------------------------------------------
Helper: validate inputs before running the main algorithm
-------------------------------------------------- */
function validateInputs(
	{
		seq,
		fwdLen,
		revLen,
		snvIndex,
		snvBase,
		tmStr,
		magnesiumStr,
		monovalentStr,
	},
	diagnostics,
) {
	const vAmp = validateAmpliconSeq(seq);
	if (!vAmp.ok) {
		diagnostics.warn('validation_failed', { check: 'amplicon_sequence' });
		goBack(vAmp.msg);
		return null;
	}

	const vPrim = validatePrimerLengths(seq.length, fwdLen, revLen);
	if (!vPrim.ok) {
		diagnostics.warn('validation_failed', { check: 'primer_lengths' });
		goBack(vPrim.msg);
		return null;
	}

	const vSnv = validateSnv(seq, fwdLen, revLen, snvIndex, snvBase);
	if (!vSnv.ok) {
		diagnostics.warn('validation_failed', { check: 'snv_selection' });
		alert(vSnv.msg);
		return null;
	}

	const vTm = validateDesiredTm(tmStr);
	if (!vTm.ok) {
		diagnostics.warn('validation_failed', { check: 'desired_tm' });
		alert(vTm.msg);
		return null;
	}

	const vConditions = validateTmConditions(magnesiumStr, monovalentStr);
	if (!vConditions.ok) {
		diagnostics.warn('validation_failed', { check: 'tm_conditions' });
		alert(vConditions.msg);
		return null;
	}

	return {
		tmC: vTm.data.tm,
		tmConditions: vConditions.data,
	};
}

function runRenderStep(diagnostics, step, render) {
	const startedAt = diagnostics.elapsedMs();
	diagnostics.info('render_step_started', { step });
	try {
		render();
		diagnostics.info('render_step_succeeded', {
			step,
			durationMs: diagnosticDuration(diagnostics, startedAt),
		});
	} catch (error) {
		diagnostics.error('render_step_failed', {
			step,
			durationMs: diagnosticDuration(diagnostics, startedAt),
			error: diagnostics.describeError(error),
		});
		throw error;
	}
}

function diagnosticDuration(diagnostics, startedAt) {
	return Math.round(Math.max(0, diagnostics.elapsedMs() - startedAt) * 10) / 10;
}

function exposeDiagnosticReport(diagnostics) {
	try {
		globalThis.uSnapbackDiagnostics = Object.freeze({
			releaseId: diagnostics.releaseId,
			runId: diagnostics.runId,
			getReport: diagnostics.getReport,
		});
	} catch {
		// The console still receives every record if the global cannot be exposed.
	}
}

export async function initResultsPage({
	diagnostics = createResultsDiagnostics(),
} = {}) {
	let stage = 'initialization';
	let overlay = null;
	let resultBox = null;
	const pageStartedAt = diagnostics.elapsedMs();
	diagnostics.installGlobalHandlers?.();
	exposeDiagnosticReport(diagnostics);
	diagnostics.info('init_started', {
		browser: detectBrowserSummary(globalThis.navigator?.userAgent),
		documentReadyState: document.readyState,
		visibilityState: document.visibilityState,
		online: globalThis.navigator?.onLine ?? null,
		secureContext: globalThis.isSecureContext ?? null,
		features: {
			worker: typeof globalThis.Worker === 'function',
			requestAnimationFrame:
				typeof globalThis.requestAnimationFrame === 'function',
		},
	});

	try {
		stage = 'dom_setup';
		const elements = {
			prevBtn: document.getElementById('prevBtn'),
			restartBtn: document.getElementById('restartBtn'),
			resultBox: document.getElementById('resultBox'),
			overlay: document.getElementById('loadingOverlay'),
			copySnapBtn: document.getElementById('copySnapSeqBtn'),
			copyLimitBtn: document.getElementById('copyLimitSeqBtn'),
			copySnapStatus: document.getElementById('copySnapStatus'),
			copyLimitStatus: document.getElementById('copyLimitStatus'),
			snapSeq: document.getElementById('snapSeq'),
			limitSeq: document.getElementById('limitSeq'),
		};
		const missingElementIds = Object.entries(elements)
			.filter(([, element]) => !element)
			.map(([id]) => id);
		if (missingElementIds.length > 0) {
			diagnostics.error('dom_check_failed', { missingElementIds });
			throw new Error(
				`Results page is missing required elements: ${missingElementIds.join(', ')}`,
			);
		}
		({ overlay, resultBox } = elements);

		wireCopyButton(elements.copySnapBtn, elements.snapSeq, {
			statusEl: elements.copySnapStatus,
		});
		wireCopyButton(elements.copyLimitBtn, elements.limitSeq, {
			statusEl: elements.copyLimitStatus,
		});

		elements.prevBtn.addEventListener('click', () => {
			window.location.href = PREV_PAGE;
		});
		elements.restartBtn.addEventListener('click', () => {
			sessionStorage.clear();
			window.location.href = START_PAGE;
		});

		stage = 'input_validation';
		const inputs = readInputs();
		diagnostics.info('inputs_read', {
			ampliconLength: inputs.seq.length,
			hasPrimerLengths:
				Number.isInteger(inputs.fwdLen) &&
				inputs.fwdLen > 0 &&
				Number.isInteger(inputs.revLen) &&
				inputs.revLen > 0,
			hasSnvSelection:
				Number.isInteger(inputs.snvIndex) &&
				typeof inputs.snvBase === 'string',
			hasDesiredTm: inputs.tmStr !== '',
			hasTmConditions:
				inputs.magnesiumStr !== '' && inputs.monovalentStr !== '',
		});
		const validated = validateInputs(inputs, diagnostics);
		if (!validated) {
			diagnostics.warn('init_stopped', { stage: 'input_validation' });
			return;
		}
		diagnostics.info(
			'validation_succeeded',
			summarizeResultsInputs(inputs, validated),
		);

		overlay.hidden = false;
		diagnostics.info('loading_overlay_shown');
		stage = 'initial_paint';
		await waitForInitialPaint({ diagnostics });

		stage = 'calculation';
		const calculationArgs = [
			inputs.seq,
			inputs.fwdLen,
			inputs.revLen,
			{ index: inputs.snvIndex, variantBase: inputs.snvBase },
			validated.tmC,
			validated.tmConditions,
		];
		const result = await calculateSnapbackForResults(calculationArgs, {
			directCalculate: createSnapback,
			diagnostics,
		});

		stage = 'rendering';
		const renderingStartedAt = diagnostics.elapsedMs();
		diagnostics.info('rendering_started');
		const snvStemIndex =
			result.descriptiveExtendedSnapback?.snvOnThreePrimeStem
				?.indexInThreePrimeStem;

		runRenderStep(diagnostics, 'stem_diagram', () => {
			renderStemDiagram(
				result.descriptiveUnExtendedSnapbackPrimer,
				result.descriptiveExtendedSnapback,
				inputs.seq[inputs.snvIndex],
				inputs.snvBase,
				snvStemIndex,
				result.matchesWild,
			);
		});
		runRenderStep(diagnostics, 'snapback_primer', () => {
			renderSnapbackPrimer(result, inputs.fwdLen, inputs.revLen);
		});
		runRenderStep(diagnostics, 'limiting_primer', () => {
			renderLimitingPrimer(result);
		});
		runRenderStep(diagnostics, 'tail_summary', () => {
			renderTailSummary(result);
		});
		runRenderStep(diagnostics, 'tm_summary', () => {
			renderTmSummary(result);
		});
		runRenderStep(diagnostics, 'stem_loop_sizes', () => {
			renderStemAndLoopSizes(result);
		});
		runRenderStep(diagnostics, 'delta_tm_table', () => {
			renderDeltaTmTable(
				result,
				inputs.seq[inputs.snvIndex],
				inputs.snvBase,
			);
		});

		resultBox.hidden = false;
		overlay.hidden = true;
		diagnostics.info('results_visible', {
			renderDurationMs: diagnosticDuration(
				diagnostics,
				renderingStartedAt,
			),
			totalDurationMs: diagnosticDuration(diagnostics, pageStartedAt),
		});
	} catch (error) {
		if (overlay) overlay.hidden = true;
		diagnostics.error('page_failed', {
			stage,
			error:
				typeof diagnostics.describeError === 'function'
					? diagnostics.describeError(error)
					: describeDiagnosticError(error),
		});
		goBack(error?.message || 'Snapback calculation failed.');
	}
}
