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
	DEFAULT_MAGNESIUM_MM,
	DEFAULT_MONOVALENT_MM,
} from '../shared/constants.js';

import {
	validateAmpliconSeq,
	validatePrimerLengths,
	validateSnv,
	validateDesiredTm,
	validateTmConditions,
} from '../shared/validators.js';

/* ------------------------ Document element and page specific constants ------------------------ */
const seq = sessionStorage.getItem('ampliconSeqCropped') || '';
const fwdLen = +sessionStorage.getItem('forwardPrimerLen') || 0;
const revLen = +sessionStorage.getItem('reversePrimerLen') || 0;
const snvIndex = +sessionStorage.getItem('snvIndex');
const snvBase = sessionStorage.getItem('snvBase');
const input = document.getElementById('desiredTm');
const magnesiumInput = document.getElementById('magnesiumMm');
const monovalentInput = document.getElementById('monovalentMm');
const form = document.getElementById('tmForm');
const prevBtn = document.getElementById('prevBtn');
const restartBtn = document.getElementById('restartBtn');
const NEXT = 'results.html'; // destination after Tm is set
const PREV = 'variant.html';

/* restore previous value */
const saved = sessionStorage.getItem('desiredTm');
if (saved) input.value = saved;

const savedMagnesium = sessionStorage.getItem('magnesiumMm');
magnesiumInput.value =
	savedMagnesium == null ? DEFAULT_MAGNESIUM_MM.toFixed(1) : savedMagnesium;

const savedMonovalent = sessionStorage.getItem('monovalentMm');
monovalentInput.value =
	savedMonovalent == null ? String(DEFAULT_MONOVALENT_MM) : savedMonovalent;

/* keep only digits on paste / typing and update storage */
input.addEventListener('input', () => {
	input.value = input.value.replace(/[^0-9]/g, '');
	sessionStorage.setItem('desiredTm', input.value);
});

magnesiumInput.addEventListener('input', () => {
	sessionStorage.setItem('magnesiumMm', magnesiumInput.value);
});

monovalentInput.addEventListener('input', () => {
	sessionStorage.setItem('monovalentMm', monovalentInput.value);
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

	const vAmp = validateAmpliconSeq(seq);
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

	const vConditions = validateTmConditions(
		magnesiumInput.value,
		monovalentInput.value,
	);
	if (!vConditions.ok) {
		alert(vConditions.msg);
		return;
	}

	sessionStorage.setItem('desiredTm', String(vTm.data.tm));
	sessionStorage.setItem(
		'magnesiumMm',
		String(vConditions.data.magnesiumMm),
	);
	sessionStorage.setItem(
		'monovalentMm',
		String(vConditions.data.monovalentMm),
	);
	window.location.href = NEXT;
});

/* back */
prevBtn.addEventListener('click', () => {
	sessionStorage.setItem('desiredTm', input.value.trim());
	sessionStorage.setItem('magnesiumMm', magnesiumInput.value.trim());
	sessionStorage.setItem('monovalentMm', monovalentInput.value.trim());
	window.location.href = PREV;
});
