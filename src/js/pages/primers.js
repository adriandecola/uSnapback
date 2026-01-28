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
	validateAmpliconSeq,
	validatePrimerRanges,
} from '../shared/validators.js';

import { getPrimerTm } from '../../script.js';

/* ------------------------ Document element and page specific constants ------------------------ */
const seq = sessionStorage.getItem('ampliconSeq') || '';
const form = document.getElementById('primerForm');
const seqBox = document.getElementById('ampliconSelect');
const pickFwdBtn = document.getElementById('pickFwdBtn');
const pickRevBtn = document.getElementById('pickRevBtn');
const clearPrimersBtn = document.getElementById('clearPrimersBtn');
const pickHint = document.getElementById('pickHint');
const previewBox = document.getElementById('primerPreview');
const fwdOut = document.getElementById('fwdPrimer');
const revOut = document.getElementById('revPrimer');
const ampLenCell = document.getElementById('ampLenCell');
const ampGcCell = document.getElementById('ampGcCell');
const fwdLenCell = document.getElementById('fwdLenCell');
const fwdGcCell = document.getElementById('fwdGcCell');
const revLenCell = document.getElementById('revLenCell');
const revGcCell = document.getElementById('revGcCell');
const fwdTmCell = document.getElementById('fwdTmCell');
const revTmCell = document.getElementById('revTmCell');
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
let fwdTmSeq = null;
let revTmSeq = null;
let fwdTmRequestId = 0;
let revTmRequestId = 0;

// Render the selectable amplicon
renderAmplicon();
setActivePick(fwdRange ? (revRange ? 'fwd' : 'rev') : 'fwd');
renderSelection();
updatePreview();
updatePrimerTms();

// Buttons
pickFwdBtn.addEventListener('click', () => setActivePick('fwd'));

pickRevBtn.addEventListener('click', () => {
	if (!fwdRange) {
		resetHint();
		return;
	}
	setActivePick('rev');
});

clearPrimersBtn.addEventListener('click', () => {
	fwdRange = null;
	revRange = null;
	clearRange('primerFwd');
	clearRange('primerRev');
	resetPrimerTms();
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
	updatePreview();
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
	updatePrimerTms();
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

function resetHint() {
	if (activePick === 'fwd') {
		setHint('Drag to select the forward primer.', false);
		return;
	}
	if (!fwdRange) {
		setHint('Select the forward primer first.', false);
		return;
	}
	setHint('Drag to select the reverse primer.', false);
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
		resetHint();
		return;
	}

	// Forward pick
	if (activePick === 'fwd') {
		// If reverse already exists, enforce ordering + gap
		if (revRange) {
			const gap = revRange.start - r.end - 1;
			if (r.end >= revRange.start) {
				alert('Forward primer must be left of reverse primer.');
				resetHint();
				return;
			}
			if (gap < MIN_GAP_BETWEEN_PRIMERS) {
				alert(
					`Need at least ${MIN_GAP_BETWEEN_PRIMERS} bases between primers.`
				);
				resetHint();
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
		resetHint();
		return;
	}

	// Must be to the right + with gap
	const gap = r.start - fwdRange.end - 1;
	if (r.start <= fwdRange.end) {
		alert('Reverse primer cannot overlap the forward primer.');
		resetHint();
		return;
	}
	if (gap < MIN_GAP_BETWEEN_PRIMERS) {
		alert(
			`Need at least ${MIN_GAP_BETWEEN_PRIMERS} bases between primers.`
		);
		resetHint();
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
	const hasSeq = Boolean(seq);
	const activeRange =
		isDragging && dragStart != null && dragEnd != null
			? normalizeRange(dragStart, dragEnd)
			: null;
	const activeLen = activeRange ? rangeLen(activeRange) : null;
	const activeSeq = activeRange ? seq.slice(activeRange.start, activeRange.end + 1) : '';
	const activeGc = activeSeq
		? ((activeSeq.match(/[GC]/g) || []).length / activeSeq.length) * 100
		: null;
	const fwdRangeForAmp =
		activePick === 'fwd' && activeRange ? activeRange : fwdRange;
	const revRangeForAmp =
		activePick === 'rev' && activeRange ? activeRange : revRange;

	// Forward primer sequence (live if dragging forward; else committed range)
	if (hasSeq && activePick === 'fwd' && activeRange) {
		const fwd = seq.slice(activeRange.start, activeRange.end + 1);
		fwdOut.textContent = fwd;
		showPreviewBox = true;
	} else if (hasSeq && fwdRange) {
		const fwd = seq.slice(fwdRange.start, fwdRange.end + 1);
		fwdOut.textContent = fwd;
		showPreviewBox = true;
	} else {
		fwdOut.textContent = '';
	}

	// Reverse primer sequence (live if dragging reverse; else committed range)
	if (hasSeq && activePick === 'rev' && activeRange) {
		const revSite = seq.slice(activeRange.start, activeRange.end + 1);
		revOut.textContent = reverseComplement(revSite);
		showPreviewBox = true;
	} else if (hasSeq && revRange) {
		const revSite = seq.slice(revRange.start, revRange.end + 1);
		revOut.textContent = reverseComplement(revSite);
		showPreviewBox = true;
	} else {
		revOut.textContent = '';
	}

	// Forward primer row (live if dragging forward; else committed range)
	if (activePick === 'fwd' && activeRange) {
		fwdLenCell.textContent = String(activeLen);
		fwdGcCell.textContent =
			activeGc == null ? '' : `${activeGc.toFixed(1)}%`;
		showPreviewBox = true;
	} else if (hasSeq && fwdRange) {
		const fwdSeq = seq.slice(fwdRange.start, fwdRange.end + 1);
		const gcPct = ((fwdSeq.match(/[GC]/g) || []).length / fwdSeq.length) * 100;
		fwdLenCell.textContent = String(fwdSeq.length);
		fwdGcCell.textContent = `${gcPct.toFixed(1)}%`;
	} else {
		fwdLenCell.textContent = '';
		fwdGcCell.textContent = '';
	}

	// Reverse primer row (live if dragging reverse; else committed range)
	if (activePick === 'rev' && activeRange) {
		revLenCell.textContent = String(activeLen);
		revGcCell.textContent =
			activeGc == null ? '' : `${activeGc.toFixed(1)}%`;
		showPreviewBox = true;
	} else if (hasSeq && revRange) {
		const revSeq = seq.slice(revRange.start, revRange.end + 1);
		const gcPct = ((revSeq.match(/[GC]/g) || []).length / revSeq.length) * 100;
		revLenCell.textContent = String(revSeq.length);
		revGcCell.textContent = `${gcPct.toFixed(1)}%`;
	} else {
		revLenCell.textContent = '';
		revGcCell.textContent = '';
	}

	// Amplicon row (span between primers, including primers)
	if (hasSeq && fwdRangeForAmp && revRangeForAmp) {
		const start = Math.min(fwdRangeForAmp.start, revRangeForAmp.start);
		const end = Math.max(fwdRangeForAmp.end, revRangeForAmp.end);
		const ampSeq = seq.slice(start, end + 1);
		const gcPct =
			((ampSeq.match(/[GC]/g) || []).length / ampSeq.length) * 100;
		ampLenCell.textContent = String(ampSeq.length);
		ampGcCell.textContent = `${gcPct.toFixed(1)}%`;
		showPreviewBox = true;
	} else {
		ampLenCell.textContent = '';
		ampGcCell.textContent = '';
	}

	previewBox.hidden = !showPreviewBox;
}

function getFwdPrimerSeq() {
	if (!fwdRange) return null;
	return seq.slice(fwdRange.start, fwdRange.end + 1);
}

function getRevPrimerSeq() {
	if (!revRange) return null;
	const revSite = seq.slice(revRange.start, revRange.end + 1);
	return reverseComplement(revSite);
}

function formatTmValue(tm) {
	if (tm == null || Number.isNaN(tm)) return '';
	return tm.toFixed(1);
}

function setTmLoading(cell) {
	if (!cell) return;
	cell.innerHTML = '';
	const spinner = document.createElement('span');
	spinner.className = 'tm-spinner';
	spinner.setAttribute('role', 'status');
	spinner.setAttribute('aria-label', 'Calculating');
	cell.appendChild(spinner);
}

function setTmValue(cell, tm) {
	if (!cell) return;
	cell.textContent = formatTmValue(tm);
}

function clearTmCell(cell) {
	if (!cell) return;
	cell.textContent = '';
}

function resetPrimerTms() {
	fwdTmSeq = null;
	revTmSeq = null;
	fwdTmRequestId += 1;
	revTmRequestId += 1;
	clearTmCell(fwdTmCell);
	clearTmCell(revTmCell);
}

async function updatePrimerTm(which, primerSeq) {
	const cell = which === 'fwd' ? fwdTmCell : revTmCell;
	if (!primerSeq) {
		if (which === 'fwd') {
			fwdTmSeq = null;
			fwdTmRequestId += 1;
		} else {
			revTmSeq = null;
			revTmRequestId += 1;
		}
		clearTmCell(cell);
		return;
	}

	const lastSeq = which === 'fwd' ? fwdTmSeq : revTmSeq;
	if (primerSeq === lastSeq && cell.textContent.trim() !== '') return;

	if (which === 'fwd') {
		fwdTmSeq = primerSeq;
		fwdTmRequestId += 1;
	} else {
		revTmSeq = primerSeq;
		revTmRequestId += 1;
	}

	const requestId = which === 'fwd' ? fwdTmRequestId : revTmRequestId;
	setTmLoading(cell);

	try {
		const tm = await getPrimerTm(primerSeq);
		const currentSeq = which === 'fwd' ? fwdTmSeq : revTmSeq;
		const currentId = which === 'fwd' ? fwdTmRequestId : revTmRequestId;
		if (requestId !== currentId || primerSeq !== currentSeq) return;
		setTmValue(cell, tm);
	} catch (err) {
		const currentId = which === 'fwd' ? fwdTmRequestId : revTmRequestId;
		if (requestId !== currentId) return;
		console.error(err);
		clearTmCell(cell);
	}
}

function updatePrimerTms() {
	if (isDragging) return;
	updatePrimerTm('fwd', getFwdPrimerSeq());
	updatePrimerTm('rev', getRevPrimerSeq());
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
	const vAmp = validateAmpliconSeq(seq);
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
	const vCropped = validateAmpliconSeq(croppedSeq, {
		messages: {
			tooShort: 'The amplicon is too short after cropping.',
		},
	});

	if (!vCropped.ok) {
		alert(vCropped.msg);
		return;
	}

	// --- Only clear downstream inputs if the primer/crop context changed ---
	const prevCropped = sessionStorage.getItem('ampliconSeqCropped');
	const prevFwdLen = sessionStorage.getItem('forwardPrimerLen');
	const prevRevLen = sessionStorage.getItem('reversePrimerLen');

	// Treat as changed only if we previously had a saved context AND it differs
	const primerContextChanged =
		(prevCropped != null && prevCropped !== vCropped.data.seq) ||
		(prevFwdLen != null && prevFwdLen !== String(fwdLen)) ||
		(prevRevLen != null && prevRevLen !== String(revLen));

	if (primerContextChanged) {
		// Dangerous downstream state
		sessionStorage.removeItem('snvIndex');

		// Also clear base so you don't keep a base with no index
		sessionStorage.removeItem('snvBase');
	}

	/* Save cropped as the canonical amplicon for the rest of the app */
	sessionStorage.setItem('ampliconSeqCropped', vCropped.data.seq);

	/* Store primer lens so downstream code can stay unchanged */
	sessionStorage.setItem('forwardPrimerLen', String(fwdLen));
	sessionStorage.setItem('reversePrimerLen', String(revLen));

	/* Store primer data */
	const forwardPrimer = seq.slice(fwdRange.start, fwdRange.end + 1);
	const reversePrimer = seq.slice(revRange.start, revRange.end + 1);

	/* Navigate forward */
	window.location.href = nextPage;
});

prevBtn.addEventListener('click', () => {
	window.location.href = prevPage;
});
