/*
  File:             renderStemDiagram.js
  Description:      DOM rendering helpers for results.html
  Author:           Adrian deCola
  Relative Path:    uSnapback/src/js/pages/renderStemDiagram.js
*/

let stemResizeObserver = null;
let stemSnvInterval = null;

// Complementary base match and helper function
const COMP = { A: 'T', T: 'A', C: 'G', G: 'C' };

const complementBase = (b) =>
	b ? COMP[String(b).toUpperCase()] || String(b).toUpperCase() : '';
// Normalize base to uppercase string
const normBase = (b) => (b ? String(b).trim().toUpperCase() : '');
// Is it a correct DNA Base?
const isDnaBase = (b) => b === 'A' || b === 'C' || b === 'G' || b === 'T';

// Auto-size the stem diagram based on its container width
function autoSizeStemDiagram(wrapper) {
	if (!wrapper) return;
	const diagram = wrapper.querySelector('.stem-diagram');
	if (!diagram) return;

	const topRow = diagram.querySelector('.stem-row--top');
	if (!topRow) return;

	const ntCount = (topRow.textContent || '').length || 1;
	const width = wrapper.clientWidth || wrapper.offsetWidth || 320;

	// Heuristic: keep between 10–22 px and roughly fit the width
	let fontSize = width / (ntCount * 1.3);
	fontSize = Math.max(10, Math.min(22, fontSize));

	wrapper.style.setProperty('--stem-font-size', `${fontSize}px`);
}

// Position the "Wild / Variant" label underneath the SNV base
function positionStemSnvLabel(wrapper) {
	if (!wrapper) return;

	const labelEl = wrapper.querySelector('#stemSnvLabel');
	const snvSpan = wrapper.querySelector('.stem-nt--snv');

	if (!labelEl || !snvSpan) return;

	const wrapperRect = wrapper.getBoundingClientRect();
	const snvRect = snvSpan.getBoundingClientRect();

	// Coordinates relative to the wrapper
	const centerX = snvRect.left + snvRect.width / 2 - wrapperRect.left;
	const labelTop = snvRect.bottom - wrapperRect.top + 6; // 6px gap under the stem

	labelEl.style.left = `${centerX}px`;
	labelEl.style.top = `${labelTop}px`;
}

function updateStemLayout(wrapper) {
	autoSizeStemDiagram(wrapper);
	positionStemSnvLabel(wrapper);
}

// Shared ResizeObserver for the stem diagram
if (typeof ResizeObserver !== 'undefined') {
	stemResizeObserver =
		stemResizeObserver ||
		new ResizeObserver((entries) => {
			entries.forEach((entry) => updateStemLayout(entry.target));
		});
}

function clearStemSnvInterval() {
	if (stemSnvInterval) {
		clearInterval(stemSnvInterval);
		stemSnvInterval = null;
	}
}

// Build a 2-row stem with angled, non-hybridizing mismatches
// at both the 5′ (loop) and 3′ (extension-block) ends.
function renderStemDiagram(
	descriptivePrimer,
	descriptiveExtended,
	inputWildBase, // from seq[snvIndex]
	inputVariantBase, // from snvBase
	snvStemIndex,
	matchesWild, // result.matchesWild
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

	if (stemResizeObserver) {
		stemResizeObserver.observe(wrapper);
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

	// Set up interval to toggle every 2 seconds
	stemSnvInterval = setInterval(() => {
		showWild = !showWild;
		applyState();
	}, 2000);
}

export { renderStemDiagram };
