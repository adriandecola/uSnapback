/*
  File:             amplicon.js
  Description:      Page logic for amplicon.html
  Author:           Adrian deCola
  Relative Path:    uSnapback/src/js/pages/amplicon.js
  Used by:          ../../pages/amplicon.html
*/

const input = document.getElementById('ampliconInput');
const form = document.getElementById('ampliconForm');
const prev = document.getElementById('prevBtn');
const statsBox = document.getElementById('ampliconStats');
const ampLenOut = document.getElementById('ampLenOut');
const ampGcOut = document.getElementById('ampGcOut');
const nextPg = 'primers.html';
const prevPg = 'start.html';
const restartBtn = document.getElementById('restartBtn');
const countBox = document.getElementById('ampliconCount');

// IME-composition guard (place near other const/let declarations)
let isComposing = false;

/* Restore previously-saved amplicon (raw, with whitespace) */
const savedRaw = sessionStorage.getItem('sequenceRaw');
if (savedRaw) input.value = savedRaw;
// Update amplicon statistics, if there is an amplicon
updateAmpliconStats();
// Update nucleotice count
updateCharCount();

/* Helper function to update amplicon statistics */
function updateAmpliconStats() {
	const raw = input.value.trim();
	const seq = raw.replace(/\s+/g, ''); // collapse whitespace
	const n = seq.length;

	if (n === 0) {
		statsBox.hidden = true;
		ampLenOut.textContent = '';
		ampGcOut.textContent = '';
		return;
	}

	const gc = (seq.match(/[GC]/g) || []).length;
	const pct = (gc / n) * 100;

	ampLenOut.textContent = String(n);
	ampGcOut.textContent = `${pct.toFixed(1)}%`;
	statsBox.hidden = false;
}

// Helper function to update the nucleotide count
function updateCharCount() {
	const seqLen = input.value.replace(/\s+/g, '').length; // nucleotides only
	countBox.textContent = `${seqLen} / ${AMPLICON_LIMIT} nt`;
	countBox.classList.toggle('over', seqLen > AMPLICON_LIMIT);
}

// Helper function to sanitize while preserving caret/selection & scroll
function sanitizeAndPreserveCaret(el) {
	// capture pre-sanitize positions
	const raw = el.value;
	const selStart = el.selectionStart;
	const selEnd = el.selectionEnd;
	const scrollY = el.scrollTop;

	// build sanitized string + mapped selection
	let out = '';
	let newStart = 0,
		newEnd = 0;

	for (let i = 0; i < raw.length; i++) {
		const ch = raw[i];
		// keep whitespace or A/C/G/T (any case); drop everything else
		const isWhite = /\s/.test(ch);
		const isBase = /[ACGTacgt]/.test(ch);
		if (isWhite || isBase) {
			out += isBase ? ch.toUpperCase() : ch;
			if (i < selStart) newStart++;
			if (i < selEnd) newEnd++;
		}
	}

	// only touch DOM if something actually changed
	if (out !== raw) {
		el.value = out;
		el.selectionStart = newStart;
		el.selectionEnd = newEnd;
		el.scrollTop = scrollY; // restore scroll in long inputs
	}
}

/* --------------------------------------------------
Event: restart button → clear storage and go to start.html
-------------------------------------------------- */
restartBtn.addEventListener('click', () => {
	sessionStorage.clear();
	window.location.href = 'start.html';
});

// IME composition listeners
input.addEventListener('compositionstart', () => {
	isComposing = true;
});
input.addEventListener('compositionend', () => {
	isComposing = false;
	sanitizeAndPreserveCaret(input);
	// keep your existing live-updates in sync
	sessionStorage.setItem('sequenceRaw', input.value.trim());
	updateAmpliconStats();
	updateCharCount();
});

input.addEventListener('input', (e) => {
	// skip during IME composition; sanitize on compositionend instead
	if (e.isComposing || isComposing) return;

	sanitizeAndPreserveCaret(input);

	// Update the raw amplicon in session storage
	sessionStorage.setItem('sequenceRaw', input.value.trim());
	// Update amplicon statistics
	updateAmpliconStats();
	// Update nucleotide count
	updateCharCount();
});

// CHANGED: Tab inserts '\t' while preserving caret, scroll, stats, and storage
input.addEventListener('keydown', (e) => {
	if (e.key !== 'Tab') return;

	// If an IME composition is active, let Tab do its default (usually focus move)
	if (isComposing) return;

	e.preventDefault();

	const s = input.selectionStart;
	const t = input.selectionEnd;
	const scrollY = input.scrollTop;

	// Insert a literal tab at the current selection
	input.value = input.value.slice(0, s) + '\t' + input.value.slice(t);

	// Place caret right after the inserted tab and restore scroll
	input.selectionStart = input.selectionEnd = s + 1;
	input.scrollTop = scrollY;

	// Sanitize + uppercase while preserving the updated caret/scroll
	sanitizeAndPreserveCaret(input);

	// Keep everything else in sync (session, stats, count)
	sessionStorage.setItem('sequenceRaw', input.value.trim());
	updateAmpliconStats();
	updateCharCount();
});

/* 2. On submit, save sequence and navigate forward */
form.addEventListener('submit', (e) => {
	e.preventDefault();
	/*---------- Checking the amplicon ----------*/
	const raw = input.value.trim(); // preserve user formatting
	const seq = raw.replace(/\s+/g, ''); // compact for validation
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

	// Store sequence
	sessionStorage.setItem('sequence', seq);
	sessionStorage.setItem('sequenceRaw', raw);

	// Navigate to next page
	window.location.href = nextPg;
});

/* 3. Handle Previous button click */
prev.addEventListener('click', () => {
	// Store sequence
	sessionStorage.setItem('sequenceRaw', input.value.trim());
	// Navigate to previous page
	window.location.href = prevPg;
});
