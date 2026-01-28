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

// Build a 2-row stem with angled, non-hybridizing mismatches
// at both the 5′ (loop) and 3′ (extension-block) ends.
function renderStemDiagram(
	descriptivePrimer,
	descriptiveExtended,
	inputWildBase, // from seq[snvIndex]
	inputVariantBase, // from snvBase
	snvStemIndex,
	matchesWild // result.matchesWild
) {
	const wrapper = document.getElementById('stemDiagramWrapper');
	const diagram = document.getElementById('stemDiagram');
	const labelEl = document.getElementById('stemSnvLabel');

	if (!wrapper || !diagram) return;
	if (!descriptivePrimer || !descriptiveExtended) return;

	// Show only the 5′ and 3′ ends of each mismatch block
	const summarizeMismatchEnds = (seq, reverseForDisplay = false) => {
		if (!seq || typeof seq !== 'string') return '';
		const trimmed = seq.trim();
		if (!trimmed) return '';

		// If the strand is being displayed in reverse orientation,
		// reverse the displayed ends so the snippet matches left→right direction.
		const first = trimmed[0];
		const last = trimmed[trimmed.length - 1];

		if (trimmed.length <= 2) {
			return reverseForDisplay
				? trimmed.split('').reverse().join('')
				: trimmed;
		}

		return reverseForDisplay ? `${last}…${first}` : `${first}…${last}`;
	};

	const reverseString = (s) =>
		typeof s === 'string' ? s.split('').reverse().join('') : '';

	const topStemSeq = descriptivePrimer.fivePrimeStem || '';
	const bottomStemSeq = descriptiveExtended.threePrimeStem || '';

	if (!topStemSeq || !bottomStemSeq) {
		wrapper.hidden = true;
		if (labelEl) {
			labelEl.textContent = '';
		}
		clearStemSnvInterval();
		return;
	}

	// Use the overlapping part of the stem on both strands
	const len = Math.min(topStemSeq.length, bottomStemSeq.length);

	// Bottom row is displayed 5′→3′ left-to-right (no change)
	const bottom = bottomStemSeq.slice(0, len);

	// Top row should be displayed 3′→5′ left-to-right.
	// Since the input is provided as 5′→3′, reverse it for display.
	const top = reverseString(topStemSeq.slice(0, len));

	// 5′ inner-loop mismatches (loop side)
	const innerLoopTopSeq =
		descriptivePrimer.fivePrimeInnerLoopMismatches || '';
	const innerLoopBottomSeq =
		descriptiveExtended.threePrimeInnerLoopMismatches || '';

	// 3′ terminal mismatches that block extension
	const terminalTopSeq =
		descriptivePrimer.fivePrimerLimSnapExtMismatches || '';
	const terminalBottomSeq =
		descriptiveExtended.threePrimerLimSnapExtMismatches || '';

	// Top row displayed 3′→5′, so reverse mismatch snippets for display.
	const leftTopDisplay = summarizeMismatchEnds(innerLoopTopSeq, true);
	const rightTopDisplay = summarizeMismatchEnds(terminalTopSeq, true);
	// Bottom row displayed 5′→3′, so leave mismatch snippets alone.
	const leftBottomDisplay = summarizeMismatchEnds(innerLoopBottomSeq, false);
	const rightBottomDisplay = summarizeMismatchEnds(terminalBottomSeq, false);

	// Validate SNV index in stem coordinates
	const validSnvIndex =
		Number.isInteger(snvStemIndex) &&
		snvStemIndex >= 0 &&
		snvStemIndex < len
			? snvStemIndex
			: null;

	// ------------------------------
	// Resolve what letters to show at the SNV on the displayed threePrimeStem.
	// ------------------------------
	const wildIn = normBase(inputWildBase);
	const varIn = normBase(inputVariantBase);

	const stemChar =
		validSnvIndex !== null ? normBase(bottom[validSnvIndex]) : '';

	let snapbackWild = '';
	let snapbackVar = '';

	if (isDnaBase(wildIn) && isDnaBase(varIn) && isDnaBase(stemChar)) {
		// The extended structure you display corresponds to the allele that the tail matches.
		// So the base we *expect* to see at the SNV on threePrimeStem is either:
		// - matchesWild ? wildIn : varIn
		// OR its complement (depending on which strand threePrimeStem represents in this case).
		const expectedMatch = matchesWild ? wildIn : varIn;

		const expectedSame = expectedMatch;
		const expectedComp = complementBase(expectedMatch);

		// Decide whether threePrimeStem is "same letters" vs "complement letters"
		// relative to your input allele bases.
		const useComplement =
			stemChar === expectedComp
				? true
				: stemChar === expectedSame
				? false
				: // fallback: if matchesWild doesn't help (unexpected data), try both possibilities
				  stemChar === complementBase(wildIn) ||
				  stemChar === complementBase(varIn);

		snapbackWild = useComplement ? complementBase(wildIn) : wildIn;
		snapbackVar = useComplement ? complementBase(varIn) : varIn;
	}

	// Clear any old content
	diagram.textContent = '';

	// ───── Top row: snapback strand ─────
	const topRow = document.createElement('div');
	topRow.className = 'stem-row stem-row--top stem-row--tail';

	if (leftTopDisplay) {
		const topLeftMismatch = document.createElement('span');
		topLeftMismatch.className =
			'stem-mismatch-block ' +
			'stem-mismatch-block--inner-loop ' +
			'stem-mismatch-block--left ' +
			'stem-mismatch-block--top';
		topLeftMismatch.textContent = leftTopDisplay;
		topRow.appendChild(topLeftMismatch);
	}

	const topCoreRow = document.createElement('span');
	topCoreRow.className = 'stem-core-row stem-core-row--top';
	topCoreRow.textContent = top;
	topRow.appendChild(topCoreRow);

	if (rightTopDisplay) {
		const topRightMismatch = document.createElement('span');
		topRightMismatch.className =
			'stem-mismatch-block ' +
			'stem-mismatch-block--terminal ' +
			'stem-mismatch-block--right ' +
			'stem-mismatch-block--top';
		topRightMismatch.textContent = rightTopDisplay;
		topRow.appendChild(topRightMismatch);
	}

	diagram.appendChild(topRow);

	// ───── Bottom row: extended product strand ─────
	const bottomRow = document.createElement('div');
	bottomRow.className = 'stem-row stem-row--bottom';

	if (leftBottomDisplay) {
		const bottomLeftMismatch = document.createElement('span');
		bottomLeftMismatch.className =
			'stem-mismatch-block ' +
			'stem-mismatch-block--inner-loop ' +
			'stem-mismatch-block--left ' +
			'stem-mismatch-block--bottom';
		bottomLeftMismatch.textContent = leftBottomDisplay;
		bottomRow.appendChild(bottomLeftMismatch);
	}

	const bottomCoreRow = document.createElement('span');
	bottomCoreRow.className = 'stem-core-row stem-core-row--bottom';

	let snvSpan = null;

	for (let i = 0; i < len; i += 1) {
		const ntSpan = document.createElement('span');
		ntSpan.className = 'stem-nt';

		if (i === validSnvIndex) {
			ntSpan.classList.add('stem-nt--snv');
			snvSpan = ntSpan;
			// content set by toggle logic below
		} else {
			ntSpan.textContent = bottom[i] || 'N';
		}

		bottomCoreRow.appendChild(ntSpan);
	}

	bottomRow.appendChild(bottomCoreRow);

	if (rightBottomDisplay) {
		const bottomRightMismatch = document.createElement('span');
		bottomRightMismatch.className =
			'stem-mismatch-block ' +
			'stem-mismatch-block--terminal ' +
			'stem-mismatch-block--right ' +
			'stem-mismatch-block--bottom';
		bottomRightMismatch.textContent = rightBottomDisplay;
		bottomRow.appendChild(bottomRightMismatch);
	}

	diagram.appendChild(bottomRow);

	// Initial sizing + observe for future resizes
	wrapper.hidden = false;
	if (labelEl) {
		labelEl.hidden = false;
	}

	// Initial sizing + label positioning
	updateStemLayout(wrapper);

	if (window._stemResizeObserver) {
		window._stemResizeObserver.observe(wrapper);
	}

	// If we don't have a valid SNV or bases, just show static stem.
	clearStemSnvInterval();

	// ---------------
	// Toggling logic of SNV
	// ---------------
	if (
		!snvSpan ||
		!labelEl ||
		!snapbackWild ||
		!snapbackVar ||
		validSnvIndex === null
	) {
		if (snvSpan && validSnvIndex !== null) {
			// If we can't resolve toggle letters, show whatever is in the stem string
			snvSpan.textContent = bottom[validSnvIndex] || 'N';
		}
		if (labelEl) labelEl.textContent = '';
		return;
	}

	const wildChar = snapbackWild;
	const variantChar = snapbackVar;

	// Start on whatever the stem actually contains (best UX)
	const stemCharNow = normBase(bottom[validSnvIndex]);
	let showWild = stemCharNow === wildChar;
	if (stemCharNow !== wildChar && stemCharNow !== variantChar) {
		// fallback if stemChar isn't a clean A/C/G/T
		showWild = true;
	}

	const applyState = () => {
		const currentBase = showWild ? wildChar : variantChar;
		const labelText = showWild
			? `Wild: ${wildChar}`
			: `Variant: ${variantChar}`;

		snvSpan.textContent = currentBase;

		snvSpan.classList.remove('stem-snv--wild', 'stem-snv--variant');
		labelEl.classList.remove('stem-snv--wild', 'stem-snv--variant');

		const cls = showWild ? 'stem-snv--wild' : 'stem-snv--variant';
		snvSpan.classList.add(cls);
		labelEl.classList.add(cls);

		labelEl.textContent = labelText;

		// Re-position label under the SNV after layout
		if (typeof requestAnimationFrame === 'function') {
			requestAnimationFrame(() => positionStemSnvLabel(wrapper));
		} else {
			positionStemSnvLabel(wrapper);
		}
	};

	// Initial state = wild
	applyState();

	window._stemSnvInterval = setInterval(() => {
		showWild = !showWild;
		applyState();
	}, 2000);
}

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
		window.location.href = 'start.html';
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
				'\n\nYou will now be redirected to the amplicon page so you can adjust your inputs.'
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

	try {
		overlay.hidden = false; // Show loading screen
		const result = await createSnapback(
			/* now referenced from global scope */
			seq,
			fwdLen,
			revLen,
			{ index: snvIndex, variantBase: snvBase },
			tmC
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
			result.matchesWild
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
