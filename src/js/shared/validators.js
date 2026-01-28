/*
  File:             validators.js
  Description:      Shared validation helpers for uSnapback page flows
  Author:           Adrian deCola
  Relative Path:    uSnapback/src/js/shared/validators.js

  Design:
  - Pure functions (no DOM, no alerts)
  - Return { ok, msg, data } so callers decide alert/redirect behavior
*/

import {
	AMPLICON_LIMIT,
	MIN_AMP_LEN,
	MIN_PRIMER_LEN,
	MIN_GAP_BETWEEN_PRIMERS,
	SNV_GAP,
	TM_MIN,
	TM_MAX,
	BASES,
} from './constants.js';

/* =========================================================
   Result helpers
========================================================= */
function ok(data = {}) {
	return { ok: true, msg: '', data };
}

function fail(msg, data = {}) {
	return { ok: false, msg, data };
}

function isInt(n) {
	return Number.isInteger(n);
}

/* Normalize a DNA sequence by removing whitespace and uppercasing */
function normalizeSeq(s) {
	if (typeof s !== 'string') return '';
	return s.replace(/\s+/g, '').toUpperCase();
}

/* =========================================================
   1) Amplicon validation
========================================================= */
export function validateAmpliconRaw(raw, opts = {}) {
	const messages = {
		missing: 'Amplicon sequence not found. Try again or restart.',
		invalidChars:
			'Amplicon contains invalid characters. Use only A, C, G, T and whitespace.',
		...opts.messages,
	};

	const r = typeof raw === 'string' ? raw.trim() : '';
	if (!r) return fail(messages.missing);

	// raw may include whitespace, but ONLY whitespace + A/C/G/T
	if (!/^[\sACGT]*$/i.test(r)) return fail(messages.invalidChars);

	const seq = r.replace(/\s+/g, '').toUpperCase();
	if (!seq) return fail(messages.missing);

	return ok({ raw: r, seq, length: seq.length });
}

export function validateAmpliconSeq(seq, opts = {}) {
	const messages = {
		missing: 'Amplicon sequence not found. Try again or restart.',
		invalidChars:
			'Amplicon must contain only A, C, G, and T (no whitespace).',
		tooShort: 'The amplicon is too short.',
		tooLong: null,
		...opts.messages,
	};

	const s = typeof seq === 'string' ? seq.trim().toUpperCase() : '';
	if (!s) return fail(messages.missing);

	// no whitespace + only A/C/G/T
	if (/\s/.test(s) || !/^[ACGT]+$/.test(s))
		return fail(messages.invalidChars);

	if (s.length < MIN_AMP_LEN) return fail(messages.tooShort);

	if (s.length > AMPLICON_LIMIT) {
		const msg =
			messages.tooLong ||
			`Amplicon exceeds ${AMPLICON_LIMIT} nucleotides (${s.length}). Please shorten it.`;
		return fail(msg, { length: s.length });
	}

	return ok({ seq: s, length: s.length });
}

/* =========================================================
   2) Primer validation by LENGTHS (post-primers page)
   Checks kept from your pages + adds the "min gap between primers"
   in a way that still works AFTER cropping:
     interiorLen = seqLen - fwdLen - revLen
========================================================= */
export function validatePrimerLengths(seqLen, fwdLen, revLen) {
	if (!isInt(seqLen) || seqLen <= 0) {
		return fail('Amplicon length is invalid.');
	}

	if (!isInt(fwdLen) || !isInt(revLen)) {
		return fail('Please enter valid primer lengths.');
	}

	if (fwdLen < MIN_PRIMER_LEN || revLen < MIN_PRIMER_LEN) {
		return fail(
			`Each primer must be at least ${MIN_PRIMER_LEN} nucleotides long.`
		);
	}

	if (fwdLen > seqLen) {
		return fail('Forward primer length exceeds amplicon length.');
	}
	if (revLen > seqLen) {
		return fail('Reverse primer length exceeds amplicon length.');
	}

	if (fwdLen + revLen > seqLen) {
		return fail('Combined primer lengths exceed amplicon length.');
	}

	/*
	  After cropping (your current flow), the primers are the ends.
	  The number of bases between them is:
	    interiorLen = seqLen - fwdLen - revLen
	  Enforce your min-gap rule consistently across pages.
	*/
	const interiorLen = seqLen - fwdLen - revLen;

	if (interiorLen < MIN_GAP_BETWEEN_PRIMERS) {
		return fail(
			`Need at least ${MIN_GAP_BETWEEN_PRIMERS} bases between primers.`
		);
	}

	return ok({
		fwdLen,
		revLen,
		interiorLen,
	});
}

/* =========================================================
   3) Primer validation by RANGES (primers.html selection)
   - both ranges exist
   - each len >= MIN_PRIMER_LEN
   - forward left of reverse
   - gap >= MIN_GAP_BETWEEN_PRIMERS
========================================================= */
function normalizeRange(r) {
	const start = Math.min(r.start, r.end);
	const end = Math.max(r.start, r.end);
	return { start, end };
}

function rangeLen(r) {
	return r.end - r.start + 1;
}

export function validatePrimerRanges(seqLen, fwdRange, revRange) {
	if (!fwdRange || !revRange) {
		return fail('Please select both forward and reverse primers.');
	}

	if (!isInt(seqLen) || seqLen <= 0) {
		return fail('Amplicon length is invalid.');
	}

	const f = normalizeRange(fwdRange);
	const r = normalizeRange(revRange);

	if (!isInt(f.start) || !isInt(f.end) || !isInt(r.start) || !isInt(r.end)) {
		return fail('Please select valid primer ranges.');
	}

	if (f.start < 0 || f.end >= seqLen || r.start < 0 || r.end >= seqLen) {
		return fail('Primer range is outside the amplicon.');
	}

	const fLen = rangeLen(f);
	const rLen = rangeLen(r);

	if (fLen < MIN_PRIMER_LEN || rLen < MIN_PRIMER_LEN) {
		return fail(
			`Each primer must be at least ${MIN_PRIMER_LEN} nucleotides long.`
		);
	}

	if (r.start <= f.end) {
		return fail('Primers cannot overlap.');
	}

	const gap = r.start - f.end - 1;

	if (gap < MIN_GAP_BETWEEN_PRIMERS) {
		return fail(
			`Need at least ${MIN_GAP_BETWEEN_PRIMERS} bases between primers.`
		);
	}

	return ok({
		fwdRange: f,
		revRange: r,
		fwdLen: fLen,
		revLen: rLen,
		gap,
	});
}

/* =========================================================
   4) SNV validation
   - index present + in range
   - base selected and valid
   - index not within primer-binding regions
   - index not too close (SNV_GAP)
   - variant base differs from wild base
========================================================= */
export function validateSnv(seq, fwdLen, revLen, snvIndex, snvBase) {
	const s = normalizeSeq(seq);

	if (!s) return fail('Amplicon sequence not found. Please restart.');

	// allow snvIndex to come in as string or number
	const rawIdx = typeof snvIndex === 'string' ? snvIndex.trim() : snvIndex;

	if (rawIdx === '') return fail('Please select a variant position.');

	const idx = typeof rawIdx === 'number' ? rawIdx : Number(rawIdx);

	if (!isInt(idx) || idx < 0 || idx >= s.length) {
		return fail('Variant index is out of range.');
	}

	if (!BASES.includes(snvBase)) {
		return fail('Please select a variant base.');
	}

	// inside primer regions?
	if (idx < fwdLen || idx >= s.length - revLen) {
		return fail(
			'SNV overlaps a primer-binding site. Choose a different index.'
		);
	}

	// too close to primers (need >= SNV_GAP bases between)
	if (idx < fwdLen + SNV_GAP) {
		return fail(
			`SNV is too close to the forward primer (need ≥ ${SNV_GAP}-bp gap).`
		);
	}

	const maxIdx = s.length - revLen - 1 - SNV_GAP;
	if (idx > maxIdx) {
		return fail(
			`SNV is too close to the reverse primer (need ≥ ${SNV_GAP}-bp gap).`
		);
	}

	// variant base must differ from wild base
	const wild = s[idx];
	if (snvBase === wild) {
		return fail(
			`Variant base must be different from the wild-type base at index ${idx} (${wild}).`
		);
	}

	return ok({ idx, wildBase: wild, variantBase: snvBase });
}

/* =========================================================
   5) Desired Tm validation
   - value present
   - whole number
   - within [TM_MIN, TM_MAX]
========================================================= */
export function validateDesiredTm(tmRaw) {
	const raw = typeof tmRaw === 'string' ? tmRaw.trim() : tmRaw;

	if (raw === '' || raw == null) return fail('Please enter a Tm value.');

	const tm = typeof raw === 'number' ? raw : Number(raw);

	if (!isInt(tm)) return fail('Tm must be a whole number.');

	if (tm < TM_MIN || tm > TM_MAX) {
		return fail(`Tm must be between ${TM_MIN} °C and ${TM_MAX} °C.`);
	}

	return ok({ tm });
}
