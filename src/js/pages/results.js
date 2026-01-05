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

// Importing validator functions
import {
	validateAmpliconSeq,
	validatePrimerLengths,
	validateSnv,
	validateDesiredTm,
} from '../shared/validators.js';

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

		/* populate UI */
		document.getElementById('snapSeq').textContent = result.snapbackSeq;
		document.getElementById('limitSeq').textContent =
			result.limitingPrimerSeq;

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
