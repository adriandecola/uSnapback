/*
  File:             resultsRender.js
  Description:      DOM rendering helpers for results.html
  Author:           Adrian deCola
  Relative Path:    uSnapback/src/js/pages/resultsRender.js
*/

function segLen(s) {
	return typeof s === 'string' ? s.trim().length : null;
}

function complementBase(base) {
	return { A: 'T', T: 'A', C: 'G', G: 'C' }[
		String(base || '').toUpperCase()
	];
}

export function renderSnapbackPrimer(result, fwdLen, revLen) {
	const snapEl = document.getElementById('snapSeq');
	const primerLabelEl = document.getElementById('snapPrimerLabel');

	// Label should reflect which primer receives the tail
	if (primerLabelEl) {
		primerLabelEl.textContent = result.tailOnForwardPrimer
			? 'Forward primer'
			: 'Reverse primer';
	}

	// The unextended primer descriptor follows the actual 5′→3′ sequence order:
	// terminal mismatch, hybridizing stem, optional loop mismatch, primer.
	const d0 = result.descriptiveUnExtendedSnapbackPrimer || {};
	const terminalMismatch = d0.fivePrimerLimSnapExtMismatches || '';
	const stemPart = d0.fivePrimeStem || '';
	const innerLoopMismatch = d0.fivePrimeInnerLoopMismatches || '';
	const describedPrimer = d0.forwardPrimer || '';
	const naturalLoopPrimerBase = describedPrimer[0] || '';
	const naturalLoopAmpliconBase =
		result.descriptiveExtendedSnapback?.stuffBetween?.slice(-1) || '';
	const hasNaturalInnerLoopMismatch = Boolean(
		!innerLoopMismatch &&
			naturalLoopPrimerBase &&
			naturalLoopAmpliconBase &&
			complementBase(naturalLoopPrimerBase) !==
				naturalLoopAmpliconBase.toUpperCase(),
	);
	const describedSnapSeq =
		terminalMismatch + stemPart + innerLoopMismatch + describedPrimer;

	const snapSeq = String(result.snapbackSeq || '');

	// Determine primer length from the chosen orientation (fallback splitting)
	let primerLen = result.tailOnForwardPrimer ? fwdLen : revLen;
	if (!Number.isInteger(primerLen) || primerLen < 1) primerLen = null;

	if (snapEl && describedSnapSeq && describedSnapSeq === snapSeq) {
		snapEl.textContent = '';

		const appendSegment = (text, className, label) => {
			if (!text) return;
			const span = document.createElement('span');
			span.className = className;
			span.textContent = text;
			span.title = label;
			snapEl.appendChild(span);
		};

		appendSegment(
			terminalMismatch,
			'seq-seg seq-seg--tail seq-seg--mismatch seq-seg--terminal-mismatch',
			'End mismatch',
		);
		appendSegment(
			stemPart,
			'seq-seg seq-seg--tail seq-seg--stem',
			'Snapback hybridizing sequence',
		);
		appendSegment(
			innerLoopMismatch,
			'seq-seg seq-seg--tail seq-seg--mismatch seq-seg--inner-loop-mismatch',
			'Added internal-loop mismatch',
		);
		if (hasNaturalInnerLoopMismatch) {
			appendSegment(
				naturalLoopPrimerBase,
				'seq-seg seq-seg--primer seq-seg--mismatch seq-seg--inner-loop-mismatch seq-seg--natural-mismatch',
				'Natural first-loop mismatch; no base was added',
			);
			appendSegment(
				describedPrimer.slice(1),
				'seq-seg seq-seg--primer',
				result.tailOnForwardPrimer ? 'Forward primer' : 'Reverse primer',
			);
		} else {
			appendSegment(
				describedPrimer,
				'seq-seg seq-seg--primer',
				result.tailOnForwardPrimer ? 'Forward primer' : 'Reverse primer',
			);
		}
		return;
	}

	let tailPart = '';
	let primerPart = '';

	if (primerLen && snapSeq.length >= primerLen) {
		primerPart = snapSeq.slice(-primerLen);
		tailPart = snapSeq.slice(0, snapSeq.length - primerLen);
	} else {
		// last resort: show plain text
		if (snapEl) snapEl.textContent = snapSeq;
	}

	// Render highlighted segments (background colors)
	if (snapEl && (tailPart || primerPart)) {
		snapEl.textContent = ''; // clear
		const tailSpan = document.createElement('span');
		tailSpan.className = 'seq-seg seq-seg--tail';
		tailSpan.textContent = tailPart;

		const primerSpan = document.createElement('span');
		primerSpan.className = 'seq-seg seq-seg--primer';
		primerSpan.textContent = primerPart;

		snapEl.appendChild(tailSpan);
		snapEl.appendChild(primerSpan);
	}
}

export function renderLimitingPrimer(result) {
	const el = document.getElementById('limitSeq');
	if (el) el.textContent = result.limitingPrimerSeq || '';
}

export function renderTailSummary(result) {
	const tailSideEl = document.getElementById('tailSide');
	const matchesWildEl = document.getElementById('matchesWild');

	if (tailSideEl) {
		tailSideEl.textContent = result.tailOnForwardPrimer
			? 'forward primer'
			: 'reverse primer';
	}

	if (matchesWildEl) {
		matchesWildEl.textContent = result.matchesWild
			? 'wild-type'
			: 'variant';
	}
}

export function renderTmSummary(result) {
	const santaWild = result.snapbackMeltingTms?.wildTm;
	document.getElementById('wildTm').textContent = Number.isFinite(santaWild)
		? santaWild.toFixed(1)
		: '—';

	const santaVariant = result.snapbackMeltingTms?.variantTm;
	document.getElementById('varTm').textContent = Number.isFinite(santaVariant)
		? santaVariant.toFixed(1)
		: '—';

	const empiricalWild = result.snapbackTmWittwer?.wildTm;
	const empiricalVariant = result.snapbackTmWittwer?.variantTm;
	const empiricalWildEl = document.getElementById('wittwerWildTm');
	const empiricalVariantEl = document.getElementById('wittwerVarTm');
	if (empiricalWildEl) {
		empiricalWildEl.textContent = Number.isFinite(empiricalWild)
			? empiricalWild.toFixed(1)
			: '—';
	}
	if (empiricalVariantEl) {
		empiricalVariantEl.textContent = Number.isFinite(empiricalVariant)
			? empiricalVariant.toFixed(1)
			: '—';
	}
}

export function renderTmConditions(tmConditions) {
	const freeMagnesiumEl = document.getElementById('freeMagnesiumMm');
	const totalMonovalentEl = document.getElementById('totalMonovalentMm');
	const freeMagnesium = Number(tmConditions?.magnesiumMm);
	const totalMonovalent = Number(tmConditions?.monovalentMm);

	if (freeMagnesiumEl) {
		freeMagnesiumEl.textContent = Number.isFinite(freeMagnesium)
			? String(freeMagnesium)
			: '—';
	}
	if (totalMonovalentEl) {
		totalMonovalentEl.textContent = Number.isFinite(totalMonovalent)
			? String(totalMonovalent)
			: '—';
	}
}

export function renderStemAndLoopSizes(result) {
	// ---------------- Stem + Loop sizes ----------------
	// Stem size (bases) = threePrimeStem length (canonical stem interval)
	const stemBases = segLen(
		result.descriptiveExtendedSnapback?.threePrimeStem
	);
	document.getElementById('stemBases').textContent = Number.isInteger(
		stemBases
	)
		? String(stemBases)
		: '—';

	// Loop size (bases) = fivePrimeInnerLoopMismatches + stuffBetween + threePrimeInnerLoopMismatches
	const fivePrimeInnerLoopMismatchesLen = segLen(
		result.descriptiveExtendedSnapback?.fivePrimeInnerLoopMismatches
	);
	const stuffBetweenLen = segLen(
		result.descriptiveExtendedSnapback?.stuffBetween
	);
	const threePrimeInnerLoopMismatchesLen = segLen(
		result.descriptiveExtendedSnapback?.threePrimeInnerLoopMismatches
	);

	const loopLength = document.getElementById('loopBases');

	// If any segment is missing, show a clear problem
	const missing = [];
	if (!Number.isInteger(fivePrimeInnerLoopMismatchesLen))
		missing.push('fivePrimeInnerLoopMismatches');
	if (!Number.isInteger(stuffBetweenLen)) missing.push('stuffBetween');
	if (!Number.isInteger(threePrimeInnerLoopMismatchesLen))
		missing.push('threePrimeInnerLoopMismatches');

	if (missing.length) {
		const msg = `Loop size error: missing/invalid segment(s): ${missing.join(
			', '
		)}.`;
		if (loopLength) {
			loopLength.textContent = '—';
			loopLength.title = msg;
		}
		console.error(msg, {
			fivePrimeInnerLoopMismatchesLen,
			stuffBetweenLen,
			threePrimeInnerLoopMismatchesLen,
		});
	} else {
		// Empty inner-loop mismatch strings are valid when no loop-side
		// mismatch is needed. Count whatever segments are present.
		const loopBases =
			fivePrimeInnerLoopMismatchesLen +
			stuffBetweenLen +
			threePrimeInnerLoopMismatchesLen;

		if (loopLength) {
			loopLength.textContent = String(loopBases);
			loopLength.title = '';
		}
	}
}

export function renderDeltaTmTable(result, wildBase, variantBase) {
	const complement = { A: 'T', T: 'A', C: 'G', G: 'C' };
	const wild = String(wildBase || '').trim().toUpperCase();
	const variant = String(variantBase || '').trim().toUpperCase();
	const setText = (id, text) => {
		const el = document.getElementById(id);
		if (el) el.textContent = text;
	};

	setText('dt-wild-heading', 'Wild-type match');
	setText('dt-var-heading', 'Variant match');
	setText('dt-fwd-heading', 'Tail on forward primer');
	setText('dt-rev-heading', 'Tail on reverse primer');

	// Populate ΔTm table (gracefully handle null/undefined)
	const fmt = (v, mismatchPair = '') => {
		if (v == null || v === '') return '—';
		const n = Number(v);
		if (!Number.isFinite(n)) return '—';
		return `${n.toFixed(1)}${mismatchPair ? ` (${mismatchPair})` : ''}`;
	};
	const d = result.meltingTempDiffs ?? {};
	const validAlleles = Boolean(complement[wild] && complement[variant]);
	const mismatchPairs = validAlleles
		? {
				forwardWild: `${variant}-${complement[wild]}`,
				forwardVariant: `${wild}-${complement[variant]}`,
				reverseWild: `${complement[variant]}-${wild}`,
				reverseVariant: `${complement[wild]}-${variant}`,
			}
		: {};

	setText(
		'dt-fwd-wild',
		fmt(d.onForwardPrimer?.matchWild, mismatchPairs.forwardWild),
	);
	setText(
		'dt-fwd-var',
		fmt(d.onForwardPrimer?.matchVariant, mismatchPairs.forwardVariant),
	);
	setText(
		'dt-rev-wild',
		fmt(d.onReversePrimer?.matchWild, mismatchPairs.reverseWild),
	);
	setText(
		'dt-rev-var',
		fmt(d.onReversePrimer?.matchVariant, mismatchPairs.reverseVariant),
	);
}
