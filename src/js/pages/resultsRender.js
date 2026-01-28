/*
  File:             resultsRender.js
  Description:      DOM rendering helpers for results.html
  Author:           Adrian deCola
  Relative Path:    uSnapback/src/js/pages/resultsRender.js
*/

function segLen(s) {
	return typeof s === 'string' ? s.trim().length : null;
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

	// Tail applies to the UNextended snapback primer: everything except the primer segment
	const d0 = result.descriptiveUnExtendedSnapbackPrimer || {};
	const expectedTail =
		(d0.fivePrimeInnerLoopMismatches || '') +
		(d0.fivePrimeStem || '') +
		(d0.fivePrimerLimSnapExtMismatches || '');

	const snapSeq = String(result.snapbackSeq || '');

	// Determine primer length from the chosen orientation (fallback splitting)
	let primerLen = result.tailOnForwardPrimer ? fwdLen : revLen;
	if (!Number.isInteger(primerLen) || primerLen < 1) primerLen = null;

	let tailPart = '';
	let primerPart = '';

	if (expectedTail && snapSeq.startsWith(expectedTail)) {
		tailPart = expectedTail;
		primerPart = snapSeq.slice(expectedTail.length);
	} else if (primerLen && snapSeq.length >= primerLen) {
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
		// If any segment length is 0, show the problem (your request)
		const zero = [];
		if (fivePrimeInnerLoopMismatchesLen === 0)
			zero.push('fivePrimeInnerLoopMismatches');
		if (stuffBetweenLen === 0) zero.push('stuffBetween');
		if (threePrimeInnerLoopMismatchesLen === 0)
			zero.push('threePrimeInnerLoopMismatches');

		if (zero.length) {
			const msg = `Loop size error: ${zero.join(
				', '
			)} length is 0 (unexpected).`;
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
}

export function renderDeltaTmTable(result) {
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
}
