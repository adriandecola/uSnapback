/*
  File:             results.js
  Description:      Page logic for results.html
  Author:           Adrian deCola
  Relative Path:    uSnapback/src/js/pages/results.js
  Used by:          ../../pages/results.html
*/
import { createSnapback } from '../script.js';

// Visualization helpers for snapback / extended products
function renderSnapbackPrimerDiagram(descriptivePrimer, descriptiveExtended) {
	const container = document.getElementById('snapSeqDiagram');
	if (!container || !descriptivePrimer) return;

	// clear any old content
	container.textContent = '';

	const segments = [
		{
			key: 'fivePrimerLimSnapExtMismatches',
			className: 'seg-mismatch-3p',
			label: '3′ blocking mismatches on tail side',
		},
		{
			key: 'fivePrimeStem',
			className: 'seg-stem',
			label: 'Stem (tail side)',
		},
		{
			key: 'fivePrimeInnerLoopMismatches',
			className: 'seg-inner-loop',
			label: 'Inner-loop mismatches on tail side',
		},
		{
			key: 'forwardPrimer',
			className: 'seg-primer',
			label: 'Forward primer',
		},
	];

	segments.forEach((seg) => {
		const seq = descriptivePrimer[seg.key];
		if (!seq || !seq.length) return;

		// one row per segment: [sequence] | [swatch + label]
		const row = document.createElement('div');
		row.className = 'sequence-row';

		const seqCell = document.createElement('div');
		seqCell.className = 'sequence-row__seq';

		const seqSpan = document.createElement('span');
		seqSpan.className = `dna-segment ${seg.className}`;
		seqSpan.textContent = seq;
		seqCell.appendChild(seqSpan);

		const labelCell = document.createElement('div');
		labelCell.className = 'sequence-row__label';

		const swatch = document.createElement('span');
		swatch.className = `legend-swatch ${seg.className}`;
		labelCell.appendChild(swatch);

		const labelText = document.createElement('span');
		labelText.textContent = seg.label; // no variable names, single line
		labelCell.appendChild(labelText);

		row.appendChild(seqCell);
		row.appendChild(labelCell);
		container.appendChild(row);
	});
}

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
	window._stemResizeObserver =
		window._stemResizeObserver ||
		new ResizeObserver((entries) => {
			entries.forEach((entry) => updateStemLayout(entry.target));
		});
}

function clearStemSnvInterval() {
	if (window._stemSnvInterval) {
		clearInterval(window._stemSnvInterval);
		window._stemSnvInterval = null;
	}
}

// Build a 2-row stem with angled, non-hybridizing mismatches
// at both the 5′ (loop) and 3′ (extension-block) ends.
function renderStemDiagram(
	descriptivePrimer,
	descriptiveExtended,
	wildBase,
	variantBase,
	snvStemIndex
) {
	const wrapper = document.getElementById('stemDiagramWrapper');
	const diagram = document.getElementById('stemDiagram');
	const labelEl = document.getElementById('stemSnvLabel');

	if (!wrapper || !diagram) return;
	if (!descriptivePrimer || !descriptiveExtended) return;

	// Show only the 5′ and 3′ ends of each mismatch block
	const summarizeMismatchEnds = (seq) => {
		if (!seq || typeof seq !== 'string') return '';
		const trimmed = seq.trim();
		if (!trimmed) return '';
		if (trimmed.length <= 2) return trimmed;
		const first = trimmed[0];
		const last = trimmed[trimmed.length - 1];
		return `${first}…${last}`;
	};

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
	const top = topStemSeq.slice(0, len);
	const bottom = bottomStemSeq.slice(0, len);

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

	const leftTopDisplay = summarizeMismatchEnds(innerLoopTopSeq);
	const leftBottomDisplay = summarizeMismatchEnds(innerLoopBottomSeq);
	const rightTopDisplay = summarizeMismatchEnds(terminalTopSeq);
	const rightBottomDisplay = summarizeMismatchEnds(terminalBottomSeq);

	// Validate SNV index in stem coordinates
	const validSnvIndex =
		Number.isInteger(snvStemIndex) &&
		snvStemIndex >= 0 &&
		snvStemIndex < len
			? snvStemIndex
			: null;

	// Clear any old content
	diagram.textContent = '';

	// ───── Top row: snapback strand ─────
	const topRow = document.createElement('div');
	topRow.className = 'stem-row stem-row--top';

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

	if (
		!snvSpan ||
		!labelEl ||
		!wildBase ||
		!variantBase ||
		validSnvIndex === null
	) {
		if (snvSpan && validSnvIndex !== null) {
			snvSpan.textContent = bottom[validSnvIndex] || 'N';
		}
		if (labelEl) {
			labelEl.textContent = '';
		}
		return;
	}

	const wildChar = wildBase;
	const variantChar = variantBase;

	let showWild = true;

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

	/* ---------- Clipboard helper ---------- */
	function copyTextToClipboard(text) {
		if (!text) return;

		// Prefer modern async clipboard API
		if (navigator.clipboard && navigator.clipboard.writeText) {
			navigator.clipboard.writeText(text).catch((err) => {
				console.error('Clipboard copy failed:', err);
			});
			return;
		}

		// Fallback for older browsers
		const textarea = document.createElement('textarea');
		textarea.value = text;
		textarea.setAttribute('readonly', '');
		textarea.style.position = 'absolute';
		textarea.style.left = '-9999px';
		document.body.appendChild(textarea);
		textarea.select();

		try {
			document.execCommand('copy');
		} catch (err) {
			console.error('execCommand copy failed:', err);
		}

		document.body.removeChild(textarea);
	}

	/* ---------- Wire up copy buttons ---------- */
	if (copySnapBtn) {
		copySnapBtn.addEventListener('click', () => {
			const el = document.getElementById('snapSeq');
			const text = el?.textContent.trim();
			copyTextToClipboard(text);
		});
	}

	if (copyLimitBtn) {
		copyLimitBtn.addEventListener('click', () => {
			const el = document.getElementById('limitSeq');
			const text = el?.textContent.trim();
			copyTextToClipboard(text);
		});
	}

	/* ---------- Nav targets ---------- */
	const PREV = 'desiredTm.html';
	const START = 'start.html';

	// maximum nt allowed
	const AMPLICON_LIMIT = 1000;

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

	const seq = sessionStorage.getItem('sequence');
	const fwdLen = +sessionStorage.getItem('forwardPrimerLen');
	const revLen = +sessionStorage.getItem('reversePrimerLen');
	const snvIndex = +sessionStorage.getItem('snvIndex');
	const snvBase = sessionStorage.getItem('snvBase');
	const tm = +sessionStorage.getItem('desiredTm');

	// Helper function
	const goBack = (msg) => {
		alert(
			msg +
				'\n\nYou will now be redirected to the amplicon page so you can adjust your inputs.'
		);
		window.location.href = 'amplicon.html';
	};

	// Type checking
	/*---------- Checking the amplicon ----------*/
	if (!seq) {
		goBack('Amplicon sequence not found. Please restart.');
		return;
	}
	if (seq.length < 33) {
		goBack('The amplicon it too short.');
		return;
	}
	if (seq.length > AMPLICON_LIMIT) {
		goBack(
			`Amplicon exceeds ${AMPLICON_LIMIT} nucleotides (${seq.length}). Please shorten it.`
		);
		return;
	}

	/*---------- Checking the primers ----------*/
	if (!Number.isInteger(fwdLen) || !Number.isInteger(revLen)) {
		goBack('Please enter valid primer lengths.');
		return;
	}
	if (fwdLen < 12 || revLen < 12) {
		goBack('Each primer must be at least 12 nucleotides long.');
		return;
	}
	if (fwdLen > seq.length) {
		// NEW
		goBack('Forward primer length exceeds amplicon length.');
		return;
	}
	if (revLen > seq.length) {
		// NEW
		goBack('Reverse primer length exceeds amplicon length.');
		return;
	}
	if (fwdLen + revLen > seq.length) {
		goBack('Combined primer lengths exceed amplicon length.');
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
	if (!Number.isInteger(tm)) {
		alert('Tm must be a whole number.');
		return;
	}
	if (tm < 40 || tm > 80) {
		alert('Tm must be between 40 °C and 80 °C.');
		return;
	}

	try {
		overlay.hidden = false; // Show loading screen
		const result = await createSnapback(
			/* now referenced from global scope */
			seq,
			fwdLen,
			revLen,
			{ index: snvIndex, variantBase: snvBase },
			tm
		);

		// For debugging
		console.log(result);

		/* populate UI */
		document.getElementById('snapSeq').textContent = result.snapbackSeq;
		document.getElementById('limitSeq').textContent =
			result.limitingPrimerSeq;
		renderSnapbackPrimerDiagram(
			result.descriptiveUnExtendedSnapbackPrimer,
			result.descriptiveExtendedSnapback
		); // segmented/color-coded parts for the snapback primer

		const snvStemIndex =
			result.descriptiveExtendedSnapback?.snvOnThreePrimeStem
				?.indexInThreePrimeStem;

		renderStemDiagram(
			result.descriptiveUnExtendedSnapbackPrimer,
			result.descriptiveExtendedSnapback,
			seq[snvIndex], // wild base from the amplicon
			snvBase, // variant base
			snvStemIndex // index within the threePrimeStem string
		); // two-row horizontal stem view with SNV toggle

		document.getElementById('tailSide').textContent =
			result.tailOnForwardPrimer ? 'Forward primer' : 'Reverse primer';
		document.getElementById('matchesWild').textContent = result.matchesWild
			? 'Wild allele'
			: 'Variant allele';

		// build display strings for Wild/Variant with labels: Wittwer, Rochester, SantaLucia
		const wittWild = result.snapbackMeltingTms?.wildTm;
		const rochWild = result.snapbackTmRochester?.wildTm;
		const slWild = result.snapbackTmSantaLucia?.wildTm;
		const wildParts = [];
		if (Number.isFinite(wittWild))
			wildParts.push(`${wittWild.toFixed(1)} -Wittwer`);
		if (Number.isFinite(rochWild))
			wildParts.push(`${rochWild.toFixed(1)} -Rochester`);
		if (Number.isFinite(slWild))
			wildParts.push(`${slWild.toFixed(1)} -SantaLucia`);
		document.getElementById('wildTm').textContent = wildParts.length
			? wildParts.join(', ')
			: '—';

		const wittVar = result.snapbackMeltingTms?.variantTm;
		const rochVar = result.snapbackTmRochester?.variantTm;
		const slVar = result.snapbackTmSantaLucia?.variantTm;
		const varParts = [];
		if (Number.isFinite(wittVar))
			varParts.push(`${wittVar.toFixed(1)} -Wittwer`);
		if (Number.isFinite(rochVar))
			varParts.push(`${rochVar.toFixed(1)} -Rochester`);
		if (Number.isFinite(slVar))
			varParts.push(`${slVar.toFixed(1)} -SantaLucia`);
		document.getElementById('varTm').textContent = varParts.length
			? varParts.join(', ')
			: '—';

		// Populate ΔTm table (gracefully handle null/undefined)
		const fmt = (v) => {
			const n = Number(v);
			return Number.isFinite(n) ? n.toFixed(1) : '—';
		};
		const d = result.meltingTempDiffs ?? {};
		document.getElementById('dt-fwd-wild').textContent = fmt(
			d.onForwardPrimer?.matchWild
		);
		document.getElementById('dt-fwd-var').textContent = fmt(
			d.onForwardPrimer?.matchVariant
		);
		document.getElementById('dt-rev-wild').textContent = fmt(
			d.onReversePrimer?.matchWild
		);
		document.getElementById('dt-rev-var').textContent = fmt(
			d.onReversePrimer?.matchVariant
		);

		/* Show results and hide loading screen */
		resultBox.hidden = false;
		overlay.hidden = true;
	} catch (err) {
		overlay.hidden = true; // Hide loading screen
		console.error(err);
		goBack(err.message || 'Snapback calculation failed.');
	}
})();
