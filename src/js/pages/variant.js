/*
  File:             variant.js
  Description:      Page logic for variant.html
  Author:           Adrian deCola
  Relative Path:    uSnapback/src/js/pages/variant.js
  Used by:          ../../pages/variant.html
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
} from '../shared/validators.js';

/* -----------------------------------------------------------------------------------------
Pull saved data from sessionStorage
+ (unary plus) converts string to a number
-------------------------------------------------------------------------------------------- */
const seq = sessionStorage.getItem('ampliconSeqCropped') || '';
const fwdLen = +sessionStorage.getItem('forwardPrimerLen') || 0;
const revLen = +sessionStorage.getItem('reversePrimerLen') || 0;

/* ------------------------ Document element and page specific constants ------------------- */
const idxIn = document.getElementById('snvIndex'); // numeric index input
const baseSel = document.getElementById('snvBase'); // <select> for variant base
const box = document.getElementById('ampliconBox'); // preview container
const form = document.getElementById('variantForm');
const prevBtn = document.getElementById('prevBtn');
const restartBtn = document.getElementById('restartBtn');
const errBox = document.getElementById('snvError');
/* Target pages */
const NEXT = 'desiredTm.html';
const PREV = 'primers.html';

/* -----------------------------------------------------------------------------------------
Restore previously entered values, if any
-------------------------------------------------------------------------------------------- */
idxIn.value = sessionStorage.getItem('snvIndex') ?? '';
baseSel.dataset.keep = sessionStorage.getItem('snvBase') ?? '';

// Use centralized setter when we have a saved index
if (idxIn.value !== '') {
	setVariantIndex(+idxIn.value); // <<< now using setter
} else {
	updateBaseOptions(); // populate <select>
	renderAmplicon(); // draw coloured sequence
}

/* --------------------------------------------------
Helper function: Index validity checker(≥ 4 bp away from either primer; i.e., 3-bp gap)
-------------------------------------------------- */
function isValidSnvIndex(i) {
	if (!Number.isInteger(i)) return false;
	if (i < 0 || i >= seq.length) return false;
	if (i < fwdLen || i >= seq.length - revLen) return false; // not inside primer regions
	if (i < fwdLen + SNV_GAP) return false; // too close to forward primer
	if (i > seq.length - revLen - 1 - SNV_GAP) return false; // too close to reverse primer
	return true;
}

/* --------------------------------------------------
Helper Function: centralized setter when choosing an index via click/keyboard
-------------------------------------------------- */
function setVariantIndex(i, focusSelect = false) {
	idxIn.value = String(i);

	// Only persist if valid under current primer/cropped context
	if (isValidSnvIndex(i)) {
		sessionStorage.setItem('snvIndex', idxIn.value);
	} else {
		sessionStorage.removeItem('snvIndex');
	}

	updateBaseOptions();
	renderAmplicon();
	if (focusSelect && !baseSel.disabled) baseSel.focus();
}

/* --------------------------------------------------
Event: restart button → clear storage and go to start.html
-------------------------------------------------- */
restartBtn.addEventListener('click', () => {
	sessionStorage.clear();
	window.location.href = 'start.html';
});

/* --------------------------------------------------
Event: user types a new index
-------------------------------------------------- */
idxIn.addEventListener('input', () => {
	/* strip any non-digit chars that slip in via paste, etc. */
	const clean = idxIn.value.replace(/[^0-9]/g, '');
	idxIn.value = clean;

	if (clean === '') {
		sessionStorage.removeItem('snvIndex');
		sessionStorage.removeItem('snvBase');
		updateBaseOptions();
		renderAmplicon();
		return;
	}
	setVariantIndex(+clean);
});

/* --------------------------------------------------
Event: click/keyboard on amplicon -> choose index
-------------------------------------------------- */
box.addEventListener('click', (e) => {
	const el = e.target.closest('span[data-idx]');
	if (!el) return;

	// don't allow choosing invalid SNV locations
	if (el.classList.contains('invalid')) return;

	setVariantIndex(+el.dataset.idx, true);
});

// optional keyboard support: Enter/Space
box.addEventListener('keydown', (e) => {
	if (
		(e.key === 'Enter' || e.key === ' ') &&
		e.target.matches('span[data-idx]')
	) {
		if (e.target.classList.contains('invalid')) return;

		e.preventDefault();
		setVariantIndex(+e.target.dataset.idx, true);
	}
});

/* --------------------------------------------------
Event: user picks a base from the <select>
-------------------------------------------------- */
baseSel.addEventListener('change', () => {
	sessionStorage.setItem('snvBase', baseSel.value || ''); // blank if none chosen
	renderAmplicon();
});

/* --------------------------------------------------
Event: back button → simply navigate
-------------------------------------------------- */
prevBtn.addEventListener('click', () => (window.location.href = PREV));

/* --------------------------------------------------
Event: Next 
-------------------------------------------------- */
form.addEventListener('submit', (e) => {
	e.preventDefault();

	const idx = +idxIn.value; // string → number
	const base = baseSel.value; // already a string

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

	const vSnv = validateSnv(seq, fwdLen, revLen, idxIn.value, baseSel.value);
	if (!vSnv.ok) {
		alert(vSnv.msg);
		return;
	}

	/* ---- Passed validation → save and continue ---- */
	sessionStorage.setItem('snvIndex', idx); // numeric ok (auto stringified)
	sessionStorage.setItem('snvBase', base);

	window.location.href = NEXT;
});

/* ==================================================
Helper: populate base <select> with the 3 alt bases
================================================== */
function updateBaseOptions() {
	const raw = idxIn.value.trim();

	/* ── 1a. no index entered yet ── */
	if (raw === '') {
		errBox.style.display = 'none'; // no error
		baseSel.style.display = 'none'; // hide selector
		baseSel.disabled = true;
		baseSel.innerHTML = '<option value="">— choose index first —</option>';
		sessionStorage.removeItem('snvBase');
		baseSel.dataset.keep = '';
		return; // done
	}

	/* ── 1b. convert and validate number ── */
	const idx = Number(raw);
	if (!Number.isInteger(idx)) return; // should never happen

	/* helper: show an error + reset selector */
	function showError(msg) {
		errBox.textContent = msg;
		errBox.style.display = 'block';
		baseSel.style.display = 'none';
		baseSel.disabled = true;
		baseSel.innerHTML = '<option value=""></option>';
		sessionStorage.removeItem('snvBase');
		baseSel.dataset.keep = '';
	}

	/* helper: build selector with 3 alt bases */
	function showSelector(refBase) {
		const opts = ['A', 'C', 'G', 'T'].filter((b) => b !== refBase);
		errBox.style.display = 'none';
		baseSel.style.display = 'inline-block';
		baseSel.disabled = false;
		baseSel.innerHTML =
			'<option value="">— choose —</option>' +
			opts.map((b) => `<option value="${b}">${b}</option>`).join('');
		if (opts.includes(baseSel.dataset.keep))
			baseSel.value = baseSel.dataset.keep;
	}

	/* range / distance checks */
	if (idx < 0 || idx >= seq.length) {
		showError('Index is outside the amplicon range.');
		return;
	}
	if (idx < fwdLen || idx >= seq.length - revLen) {
		showError('SNV overlaps a primer-binding site.');
		return;
	}
	if (idx < fwdLen + SNV_GAP || idx > seq.length - revLen - 1 - SNV_GAP) {
		showError('SNV is too close to a primer-binding site.');
		return;
	}

	/* valid position → populate selector */
	showSelector(seq[idx]);
}

/* ==================================================
Helper: render coloured amplicon
================================================== */
function renderAmplicon() {
	const rawIdx = idxIn.value.trim();
	const idx = rawIdx === '' ? -1 : +rawIdx; // -1 -> “no highlight”
	const base = baseSel.value;
	let html = '';

	for (let i = 0; i < seq.length; i++) {
		let classes = ['nt'];
		let nt = seq[i];

		// primer highlight (push correct class)
		if (i < fwdLen || i >= seq.length - revLen) classes.push('primer');

		// SNV highlight (push correct class)
		if (i === idx) {
			classes.push('snv');
			nt = base && 'ACGT'.includes(base) ? `${seq[i]}/${base}` : seq[i];
		}

		// Mark bases that are not valid clickable SNV positions
		if (!isValidSnvIndex(i)) classes.push('invalid');

		// Add data-idx + tabindex so each base is clickable/focusable and then push classes
		html += `<span data-idx="${i}" class="${classes.join(
			' '
		)}" role="button" tabindex="0" title="idx ${i}: ${
			seq[i]
		}">${nt}</span>`;
	}
	box.innerHTML = html;
}
