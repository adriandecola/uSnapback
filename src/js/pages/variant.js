/*
  File:             variant.js
  Description:      Page logic for variant.html
  Author:           Adrian deCola
  Relative Path:    uSnapback/src/js/pages/variant.js
  Used by:          ../../pages/variant.html
*/

/* --------------------------------------------------
Pull saved data from sessionStorage
+ (unary plus) converts string → number
-------------------------------------------------- */
const raw = sessionStorage.getItem('sequenceRaw') || '';
const seq = sessionStorage.getItem('sequence') || '';
const fwdLen = +sessionStorage.getItem('forwardPrimerLen') || 0;
const revLen = +sessionStorage.getItem('reversePrimerLen') || 0;

/* --------------------------------------------------
Cache DOM elements
-------------------------------------------------- */
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

/* --------------------------------------------------
Restore previously entered values, if any
-------------------------------------------------- */
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
	if (i < fwdLen + 3) return false; // too close to forward primer
	if (i > seq.length - revLen - 4) return false; // too close to reverse primer
	return true;
}

/* --------------------------------------------------
Helper Function: centralized setter when choosing an index via click/keyboard
-------------------------------------------------- */
function setVariantIndex(i, focusSelect = false) {
	idxIn.value = String(i);
	sessionStorage.setItem('snvIndex', idxIn.value);
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
Event: click/keyboard on amplicon → choose index
-------------------------------------------------- */
// NEW
box.addEventListener('click', (e) => {
	const el = e.target.closest('span[data-idx]');
	if (!el) return;
	setVariantIndex(+el.dataset.idx, true);
});

// NEW (optional keyboard support: Enter/Space)
box.addEventListener('keydown', (e) => {
	if (
		(e.key === 'Enter' || e.key === ' ') &&
		e.target.matches('span[data-idx]')
	) {
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
	if (idxIn.value.trim() === '') {
		alert('Please enter a variant index.');
		return;
	}
	if (!Number.isInteger(idx) || idx < 0 || idx >= seq.length) {
		alert('Variant index is out of range.');
		return;
	}
	if (!'ACGT'.includes(base)) {
		alert('Please select a variant base.');
		return;
	}
	/* overlaps either primer region */
	if (idx < fwdLen || idx >= seq.length - revLen) {
		alert('SNV overlaps a primer-binding site. Choose a different index.');
		return;
	}
	/* too close—need ≥ 3-bp gap (= 4 bp away) */
	if (idx < fwdLen + 3) {
		alert('SNV is too close to the forward primer (need ≥ 3-bp gap).');
		return;
	}
	if (idx > seq.length - revLen - 4) {
		alert('SNV is too close to the reverse primer (need ≥ 3-bp gap).');
		return;
	}
	// make sure variant base is different from wild-type base
	if (base === seq[idx]) {
		alert(
			`Variant base must be different from the wild-type base at index ${idx} (${seq[idx]}).`
		);
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
	if (idx < fwdLen) {
		showError('SNV overlaps a primer-binding site.');
		return;
	}
	if (idx < fwdLen + 3 || idx > seq.length - revLen - 4) {
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
