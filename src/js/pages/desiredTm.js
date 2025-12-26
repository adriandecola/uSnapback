/*
  File:             desiredTm.js
  Description:      Page logic for desiredTm.html
  Author:           Adrian deCola
  Relative Path:    uSnapback/src/js/pages/desiredTm.js
  Used by:          ../../pages/desiredTm.html
*/

/* ---------------------------------------- Imports --------------------------------------- */
import {
	AMPLICON_LIMIT,
	MIN_AMP_LEN,
	MIN_PRIMER_LEN,
	MIN_GAP_BETWEEN_PRIMERS,
	SNV_GAP,
	TM_MIN,
	TM_MAX,
	BASES,
} from '../shared/constants.js';

import {
	validateAmplicon,
	validatePrimerLengths,
	validateSnv,
	validateDesiredTm,
} from '../shared/validators.js';

/* ------------------------ Document element and page specific constants ------------------------ */
const raw = sessionStorage.getItem('sequenceRaw') || '';
const seq = sessionStorage.getItem('sequence') || '';
const fwdLen = +sessionStorage.getItem('forwardPrimerLen') || 0;
const revLen = +sessionStorage.getItem('reversePrimerLen') || 0;
const snvIndex = +sessionStorage.getItem('snvIndex');
const snvBase = sessionStorage.getItem('snvBase');
const input = document.getElementById('desiredTm');
const form = document.getElementById('tmForm');
const prevBtn = document.getElementById('prevBtn');
const restartBtn = document.getElementById('restartBtn');
const NEXT = 'results.html'; // destination after Tm is set
const PREV = 'variant.html';

/* restore previous value */
const saved = sessionStorage.getItem('desiredTm');
if (saved) input.value = saved;

/* keep only digits on paste / typing and update storage */
input.addEventListener('input', () => {
	input.value = input.value.replace(/[^0-9]/g, '');
	sessionStorage.setItem('desiredTm', input.value);
});

/* --------------------------------------------------
Event: restart button → clear storage and go to start.html
-------------------------------------------------- */
restartBtn.addEventListener('click', () => {
	sessionStorage.clear();
	window.location.href = 'start.html';
});

/* submit / next */
form.addEventListener('submit', (e) => {
	e.preventDefault();
	const tm = +input.value;

	const vAmp = validateAmplicon(seq);
	if (!vAmp.ok) {
		alert(vAmp.msg);
		return;
	}

	const vPrim = validatePrimerLengths(seq.length, fwdLen, revLen);
	if (!vPrim.ok) {
		alert(vPrim.msg);
		return;
	}

	const vSnv = validateSnv(seq, fwdLen, revLen, snvIndex, snvBase);
	if (!vSnv.ok) {
		alert(vSnv.msg);
		return;
	}

	const vTm = validateDesiredTm(input.value);
	if (!vTm.ok) {
		alert(vTm.msg);
		return;
	}

	sessionStorage.setItem('desiredTm', String(vTm.data.tm));
	window.location.href = NEXT;
});

/* back */
prevBtn.addEventListener('click', () => {
	sessionStorage.setItem('desiredTm', input.value.trim());
	window.location.href = PREV;
});
