/*
  File:             primers.js
  Description:      Page logic for primers.html
  Author:           Adrian deCola
  Relative Path:    uSnapback/src/js/pages/primers.js
  Used by:          ../../pages/primers.html
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
	validatePrimerRanges,
} from '../shared/validators.js';

/* ------------------------ Document element and page specific constants ------------------------ */
const raw = sessionStorage.getItem('sequenceRaw') || '';
const seq = sessionStorage.getItem('sequence') || '';
const form = document.getElementById('primerForm');
const seqBox = document.getElementById('ampliconSelect');
const pickFwdBtn = document.getElementById('pickFwdBtn');
const pickRevBtn = document.getElementById('pickRevBtn');
const clearPrimersBtn = document.getElementById('clearPrimersBtn');
const pickHint = document.getElementById('pickHint');
const previewBox = document.getElementById('primerPreview');
const fwdOut = document.getElementById('fwdPrimer');
const revOut = document.getElementById('revPrimer');
const fwdLine = fwdOut.parentElement;
const revLine = revOut.parentElement;
const prevBtn = document.getElementById('prevBtn');
const restartBtn = document.getElementById('restartBtn');
const nextPage = 'variant.html';
const prevPage = 'amplicon.html';

let activePick = 'fwd'; // 'fwd' | 'rev'
let isDragging = false;
let dragStart = null;
let dragEnd = null;

let fwdRange = loadRange('primerFwd', seq.length);
let revRange = loadRange('primerRev', seq.length);

// Render the selectable amplicon
renderAmplicon();
setActivePick(fwdRange ? (revRange ? 'fwd' : 'rev') : 'fwd');
renderSelection();
updatePreview();

// Buttons
pickFwdBtn.addEventListener('click', () => setActivePick('fwd'));

pickRevBtn.addEventListener('click', () => {
	if (!fwdRange) {
		setHint('Select the forward primer first.', true);
		return;
	}
	setActivePick('rev');
});

clearPrimersBtn.addEventListener('click', () => {
	fwdRange = null;
	revRange = null;
	clearRange('primerFwd');
	clearRange('primerRev');
	setActivePick('fwd');
	renderSelection();
	updatePreview();
	setHint('Cleared. Select the forward primer.', false);
});

// Pointer selection (drag)
seqBox.addEventListener('pointerdown', (e) => {
	const nt = e.target.closest('.nt');
	if (!nt) return;

	// don't allow starting a drag on an invalid base while picking reverse
	if (activePick === 'rev' && nt.classList.contains('sel-invalid')) return;

	isDragging = true;
	dragStart = Number(nt.dataset.i);
	dragEnd = dragStart;

	seqBox.setPointerCapture(e.pointerId);
	renderSelection(); // show drag highlight immediately
});

seqBox.addEventListener('pointermove', (e) => {
	if (!isDragging) return;

	const el = document.elementFromPoint(e.clientX, e.clientY);
	const nt = el && el.closest ? el.closest('.nt') : null;
	if (!nt || !seqBox.contains(nt)) return;

	dragEnd = Number(nt.dataset.i);
	renderSelection();
});

seqBox.addEventListener('pointerup', () => {
	if (!isDragging) return;
	isDragging = false;

	const range = normalizeRange(dragStart, dragEnd);
	applyPickedRange(range);
	dragStart = null;
	dragEnd = null;

	renderSelection();
	updatePreview();
});

seqBox.addEventListener('pointercancel', () => {
	isDragging = false;
	dragStart = null;
	dragEnd = null;
	renderSelection();
});

function renderAmplicon() {
	if (!seq) {
		seqBox.textContent = '';
		return;
	}
	seqBox.innerHTML = '';
	const frag = document.createDocumentFragment();
	for (let i = 0; i < seq.length; i++) {
		const span = document.createElement('span');
		span.className = 'nt';
		span.dataset.i = String(i);
		span.textContent = seq[i];
		frag.appendChild(span);
	}
	seqBox.appendChild(frag);
}

function normalizeRange(a, b) {
	const start = Math.min(a, b);
	const end = Math.max(a, b);
	return { start, end };
}

function rangeLen(r) {
	return r.end - r.start + 1;
}

function setActivePick(which) {
	activePick = which;
	pickFwdBtn.classList.toggle('mini-btn--active', which === 'fwd');
	pickRevBtn.classList.toggle('mini-btn--active', which === 'rev');

	if (which === 'fwd') {
		setHint('Drag to select the forward primer.', false);
	} else {
		setHint('Drag to select the reverse primer.', false);
	}
	renderSelection();
}

function setHint(msg, isError) {
	pickHint.textContent = msg;
	pickHint.classList.toggle('error', Boolean(isError));
}

function loadRange(prefix, seqLen) {
	const s = sessionStorage.getItem(prefix + 'Start');
	const e = sessionStorage.getItem(prefix + 'End');
	if (s == null || e == null) return null;

	const start = Number(s);
	const end = Number(e);
	if (!Number.isInteger(start) || !Number.isInteger(end)) return null;
	if (start < 0 || end < 0 || start >= seqLen || end >= seqLen) return null;

	const r = normalizeRange(start, end);
	return r;
}

function saveRange(prefix, r) {
	sessionStorage.setItem(prefix + 'Start', String(r.start));
	sessionStorage.setItem(prefix + 'End', String(r.end));
}

function clearRange(prefix) {
	sessionStorage.removeItem(prefix + 'Start');
	sessionStorage.removeItem(prefix + 'End');
}

function applyPickedRange(r) {
	const len = rangeLen(r);
	if (len < MIN_PRIMER_LEN) {
		alert(`Primer must be at least ${MIN_PRIMER_LEN} nucleotides.`);
		setHint(`Too short (${len}). Min is ${MIN_PRIMER_LEN}.`, true);
		return;
	}

	// Forward pick
	if (activePick === 'fwd') {
		// If reverse already exists, enforce ordering + gap
		if (revRange) {
			const gap = revRange.start - r.end - 1;
			if (r.end >= revRange.start) {
				alert('Forward primer must be left of reverse primer.');
				setHint('Forward must be left of reverse.', true);
				return;
			}
			if (gap < MIN_GAP_BETWEEN_PRIMERS) {
				alert(
					`Need at least ${MIN_GAP_BETWEEN_PRIMERS} bases between primers.`
				);
				setHint('Not enough space between primers.', true);
				return;
			}
		}

		fwdRange = r;
		saveRange('primerFwd', fwdRange);

		// Auto-switch to reverse selection
		setActivePick('rev');
		return;
	}

	// Reverse pick
	if (!fwdRange) {
		alert('Select the forward primer first.');
		setHint('Select the forward primer first.', true);
		return;
	}

	// Must be to the right + with gap
	const gap = r.start - fwdRange.end - 1;
	if (r.start <= fwdRange.end) {
		alert('Reverse primer cannot overlap the forward primer.');
		setHint('Reverse overlaps forward.', true);
		return;
	}
	if (gap < MIN_GAP_BETWEEN_PRIMERS) {
		alert(
			`Need at least ${MIN_GAP_BETWEEN_PRIMERS} bases between primers.`
		);
		setHint(`Gap is ${gap}. Min is ${MIN_GAP_BETWEEN_PRIMERS}.`, true);
		return;
	}

	revRange = r;
	saveRange('primerRev', revRange);
	setHint('Both primers selected. Click Next.', false);
}

function renderSelection() {
	const nts = seqBox.querySelectorAll('.nt');
	if (!nts.length) return;

	nts.forEach((n) =>
		n.classList.remove('sel-fwd', 'sel-rev', 'sel-drag', 'sel-invalid')
	);

	// Mark invalid region while picking reverse (anything too far left)
	if (activePick === 'rev' && fwdRange) {
		const invalidMax = Math.min(
			fwdRange.end + MIN_GAP_BETWEEN_PRIMERS,
			nts.length - 1
		);
		for (let i = 0; i <= invalidMax; i++)
			nts[i].classList.add('sel-invalid');
	}

	// Final selections
	if (fwdRange) {
		for (let i = fwdRange.start; i <= fwdRange.end; i++)
			nts[i].classList.add('sel-fwd');
	}
	if (revRange) {
		for (let i = revRange.start; i <= revRange.end; i++)
			nts[i].classList.add('sel-rev');
	}

	// Drag preview (overlay-ish)
	if (isDragging && dragStart != null && dragEnd != null) {
		const r = normalizeRange(dragStart, dragEnd);
		for (let i = r.start; i <= r.end; i++) nts[i].classList.add('sel-drag');
	}
}

/* --- Helper Function: update preview if applicable--- */
function updatePreview() {
	let showPreviewBox = false;

	if (seq && fwdRange) {
		const fwd = seq.slice(fwdRange.start, fwdRange.end + 1);
		fwdOut.textContent = fwd;
		fwdLine.style.display = 'block';
		showPreviewBox = true;
	} else {
		fwdOut.textContent = '';
		fwdLine.style.display = 'none';
	}

	if (seq && revRange) {
		const revSite = seq.slice(revRange.start, revRange.end + 1);
		revOut.textContent = reverseComplement(revSite);
		revLine.style.display = 'block';
		showPreviewBox = true;
	} else {
		revOut.textContent = '';
		revLine.style.display = 'none';
	}

	previewBox.hidden = !showPreviewBox;
}

/* Helper Function: return reverse-complement of sequence */
function reverseComplement(dna) {
	const comp = { A: 'T', T: 'A', C: 'G', G: 'C' };
	return [...dna]
		.reverse()
		.map((b) => comp[b] || '?')
		.join('');
}

/* --------------------------------------------------
Event: restart button → clear storage and go to start.html
-------------------------------------------------- */
restartBtn.addEventListener('click', () => {
	sessionStorage.clear();
	window.location.href = 'start.html';
});

/* --- Handle form submission (Next) --- */
form.addEventListener('submit', (e) => {
	e.preventDefault();

	// Validate uncropped amplicon
	const vAmp = validateAmplicon(seq);
	if (!vAmp.ok) {
		alert(vAmp.msg);
		return;
	}

	// Validate primers
	const vRanges = validatePrimerRanges(seq.length, fwdRange, revRange);
	if (!vRanges.ok) {
		alert(vRanges.msg);
		return;
	}

	const { fwdLen, revLen } = vRanges.data;

	/* Crop the amplicon so primers become the ends (keeps downstream logic the same) */
	const croppedSeq = seq.slice(fwdRange.start, revRange.end + 1);

	// Revalidate cropped amplicon
	const vCropped = validateAmplicon(croppedSeq, {
		messages: {
			tooShort: 'The amplicon is too short after cropping.',
		},
	});
	if (!vCropped.ok) {
		alert(vCropped.msg);
		return;
	}

	/* Preserve the original (so Back to amplicon.html shows what user typed) */
	if (!sessionStorage.getItem('sequenceFull')) {
		sessionStorage.setItem('sequenceFull', seq);
		sessionStorage.setItem('sequenceRawFull', raw || seq);
	}

	/* Save cropped as the canonical amplicon for the rest of the app */
	sessionStorage.setItem('sequence', croppedSeq);
	sessionStorage.setItem('sequenceRaw', croppedSeq);

	/* Store primer lens so downstream code can stay unchanged */
	sessionStorage.setItem('forwardPrimerLen', String(fwdLen));
	sessionStorage.setItem('reversePrimerLen', String(revLen));

	/* Rebase ranges so if user comes back here, primers show at the ends */
	saveRange('primerFwd', { start: 0, end: fwdLen - 1 });
	saveRange('primerRev', {
		start: croppedSeq.length - revLen,
		end: croppedSeq.length - 1,
	});

	/* Store primer data */
	const forwardPrimer = seq.slice(fwdRange.start, fwdRange.end + 1);
	const reversePrimer = seq.slice(revRange.start, revRange.end + 1);

	sessionStorage.setItem('forwardPrimerLen', fwdLen);
	sessionStorage.setItem('reversePrimerLen', revLen);

	/* Navigate forward */
	window.location.href = nextPage;
});

prevBtn.addEventListener('click', () => {
	// Restore full amplicon if we previously cropped it
	const full = sessionStorage.getItem('sequenceFull');
	const rawFull = sessionStorage.getItem('sequenceRawFull');

	if (full) sessionStorage.setItem('sequence', full);
	if (rawFull) sessionStorage.setItem('sequenceRaw', rawFull);

	window.location.href = prevPage;
});
