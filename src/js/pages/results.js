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

(async () => {
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

	/* ---------- Nav targets ---------- */
	const PREV = 'desiredTm.html';
	const START = 'start.html';

	/* ---------- Back button ---------- */
	prevBtn.addEventListener('click', () => {
		/* keep current inputs intact – just step back */
		window.location.href = PREV;
	});

	/* --------------------------------------------------
    Event: restart button → clear storage and go to start.html
    -------------------------------------------------- */
	restartBtn.addEventListener('click', () => {
		sessionStorage.clear();
		window.location.href = START;
	});

	/* ---------- Pull inputs ------------ */

	const seq = sessionStorage.getItem('ampliconSeqCropped') || '';
	const fwdLen = +sessionStorage.getItem('forwardPrimerLen');
	const revLen = +sessionStorage.getItem('reversePrimerLen');
	const snvIndex = +sessionStorage.getItem('snvIndex');
	const snvBase = sessionStorage.getItem('snvBase');
	const tmStr = sessionStorage.getItem('desiredTm') ?? '';

	// Helper function
	const goBack = (msg) => {
		alert(
			msg +
				'\n\nYou will now be redirected to the amplicon page so you can adjust your inputs.',
		);
		window.location.href = 'amplicon.html';
	};

	// Validating Inputs
	const vAmp = validateAmpliconSeq(seq);
	if (!vAmp.ok) {
		goBack(vAmp.msg);
		return;
	}

	const vPrim = validatePrimerLengths(seq.length, fwdLen, revLen);
	if (!vPrim.ok) {
		goBack(vPrim.msg);
		return;
	}

	const vSnv = validateSnv(seq, fwdLen, revLen, snvIndex, snvBase);
	if (!vSnv.ok) {
		alert(vSnv.msg);
		return;
	}

	const vTm = validateDesiredTm(tmStr);
	if (!vTm.ok) {
		alert(vTm.msg);
		return;
	}
	const tmC = vTm.data.tm;

	/* ---------- Main logic ---------- */
	try {
		overlay.hidden = false; // Show loading screen
		const result = await createSnapback(
			/* now referenced from global scope */
			seq,
			fwdLen,
			revLen,
			{ index: snvIndex, variantBase: snvBase },
			tmC,
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
			seq[snvIndex],
			snvBase,
			snvStemIndex,
			result.matchesWild,
		);

		renderSnapbackPrimer(result, fwdLen, revLen);
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
})();
