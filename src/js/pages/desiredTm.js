/*
  File:             desiredTm.js
  Description:      Page logic for desiredTm.html
  Author:           Adrian deCola
  Relative Path:    uSnapback/src/js/pages/desiredTm.js
  Used by:          ../../pages/desiredTm.html
*/

import {
	AMPLICON_LIMIT,
	MIN_AMP_LEN,
	MIN_PRIMER_LEN,
	MIN_GAP_BETWEEN_PRIMERS,
	SNV_GAP,
	TM_MIN,
	TM_MAX,
} from '../shared/constants.js';

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

	/*---------- Checking the amplicon ----------*/
	if (!seq) {
		alert('Amplicon sequence not found. Please restart.');
		return;
	}
	if (seq.length < 33) {
		alert('The amplicon it too short.');
		return;
	}
	if (seq.length > AMPLICON_LIMIT) {
		alert(
			`Amplicon exceeds ${AMPLICON_LIMIT} nucleotides (${seq.length}). Please shorten it.`
		);
		return;
	}

	/*---------- Checking the primers ----------*/
	if (!Number.isInteger(fwdLen) || !Number.isInteger(revLen)) {
		alert('Please enter valid primer lengths.');
		return;
	}
	if (fwdLen < 12 || revLen < 12) {
		alert('Each primer must be at least 12 nucleotides long.');
		return;
	}
	if (fwdLen > seq.length) {
		// NEW
		alert('Forward primer length exceeds amplicon length.');
		return;
	}
	if (revLen > seq.length) {
		// NEW
		alert('Reverse primer length exceeds amplicon length.');
		return;
	}
	if (fwdLen + revLen > seq.length) {
		alert('Combined primer lengths exceed amplicon length.');
		return;
	}

	/*---------- Checking variant index and base ----------*/
	if (!Number.isInteger(snvIndex) || snvIndex < 0 || snvIndex >= seq.length) {
		alert('Variant index is out of range.');
		return;
	}
	if (!'ACGT'.includes(snvBase)) {
		alert('Please select a variant base.');
		return;
	}
	/* overlaps either primer region */
	if (snvIndex < fwdLen || snvIndex >= seq.length - revLen) {
		alert('SNV overlaps a primer-binding site. Choose a different index.');
		return;
	}
	/* too close—need ≥ 3-bp gap (= 4 bp away) */
	if (snvIndex < fwdLen + 3) {
		alert('SNV is too close to the forward primer (need ≥ 3-bp gap).');
		return;
	}
	if (snvIndex > seq.length - revLen - 4) {
		alert('SNV is too close to the reverse primer (need ≥ 3-bp gap).');
		return;
	}
	// make sure variant base is different from wild-type base
	if (snvBase === seq[snvIndex]) {
		alert(
			`Variant base must be different from the wild-type base at index ${snvIndex} (${seq[snvIndex]}).`
		);
		return;
	}

	/*---------- Checking Tm ----------*/
	if (input.value.trim() === '') {
		alert('Please enter a Tm value.');
		return;
	}
	if (!Number.isInteger(tm)) {
		alert('Tm must be a whole number.');
		return;
	}
	if (tm < 40 || tm > 80) {
		alert('Tm must be between 40 °C and 80 °C.');
		return;
	}

	sessionStorage.setItem('desiredTm', tm);
	window.location.href = NEXT;
});

/* back */
prevBtn.addEventListener('click', () => {
	sessionStorage.setItem('desiredTm', input.value.trim());
	window.location.href = PREV;
});
