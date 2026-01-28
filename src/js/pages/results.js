/*
  File:             results.js
  Description:      Page logic for results.html
  Author:           Adrian deCola
  Relative Path:    uSnapback/src/js/pages/results.js
  Used by:          ../../pages/results.html
*/

/* ---------------------------------------- Imports --------------------------------------- */
import { createSnapback } from '../../script.js';
// Dont import constants here as all inputs have already been collected
import { wireCopyButton } from '../shared/clipboard.js';
import { renderStemDiagram } from './resultsStemDiagram.js';
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
	};
}

/* --------------------------------------------------
Helper: validate inputs before running the main algorithm
-------------------------------------------------- */
function validateInputs({ seq, fwdLen, revLen, snvIndex, snvBase, tmStr }) {
	const vAmp = validateAmpliconSeq(seq);
	if (!vAmp.ok) {
		goBack(vAmp.msg);
		return null;
	}

	const vPrim = validatePrimerLengths(seq.length, fwdLen, revLen);
	if (!vPrim.ok) {
		goBack(vPrim.msg);
		return null;
	}

	const vSnv = validateSnv(seq, fwdLen, revLen, snvIndex, snvBase);
	if (!vSnv.ok) {
		alert(vSnv.msg);
		return null;
	}

	const vTm = validateDesiredTm(tmStr);
	if (!vTm.ok) {
		alert(vTm.msg);
		return null;
	}

	return { tmC: vTm.data.tm };
}

async function initResultsPage() {
	/* ---------- DOM elements ---------- */
	const prevBtn = document.getElementById('prevBtn');
	const restartBtn = document.getElementById('restartBtn');
	const resultBox = document.getElementById('resultBox');
	const overlay = document.getElementById('loadingOverlay');
	/* Copy buttons */
	const copySnapBtn = document.getElementById('copySnapSeqBtn');
	const copyLimitBtn = document.getElementById('copyLimitSeqBtn');

	/* ---------- Wire up copy buttons ---------- */
	wireCopyButton(copySnapBtn, document.getElementById('snapSeq'));
	wireCopyButton(copyLimitBtn, document.getElementById('limitSeq'));

	/* ---------- Back button ---------- */
	prevBtn.addEventListener('click', () => {
		/* keep current inputs intact – just step back */
		window.location.href = PREV_PAGE;
	});

	/* --------------------------------------------------
	Event: restart button → clear storage and go to start.html
	-------------------------------------------------- */
	restartBtn.addEventListener('click', () => {
		sessionStorage.clear();
		window.location.href = START_PAGE;
	});

	/* ---------- Pull and validate inputs ---------- */
	const inputs = readInputs();
	const validated = validateInputs(inputs);
	if (!validated) return;

	/* ---------- Main logic ---------- */
	try {
		overlay.hidden = false; // Show loading screen during compute/render
		const result = await createSnapback(
			inputs.seq,
			inputs.fwdLen,
			inputs.revLen,
			{ index: inputs.snvIndex, variantBase: inputs.snvBase },
			validated.tmC,
		);

		// For debugging
		console.log(result);

		const snvStemIndex =
			result.descriptiveExtendedSnapback?.snvOnThreePrimeStem
				?.indexInThreePrimeStem;

		// two-row horizontal stem view with SNV toggle
		renderStemDiagram(
			result.descriptiveUnExtendedSnapbackPrimer,
			result.descriptiveExtendedSnapback,
			inputs.seq[inputs.snvIndex],
			inputs.snvBase,
			snvStemIndex,
			result.matchesWild,
		);

		renderSnapbackPrimer(result, inputs.fwdLen, inputs.revLen);
		renderLimitingPrimer(result);

		renderTailSummary(result);
		renderTmSummary(result);
		renderStemAndLoopSizes(result);
		renderDeltaTmTable(result);

		/* Show results and hide loading screen */
		resultBox.hidden = false;
		overlay.hidden = true;
	} catch (err) {
		overlay.hidden = true; // Hide loading screen
		console.error(err);
		goBack(err.message || 'Snapback calculation failed.');
	}
}

initResultsPage();
