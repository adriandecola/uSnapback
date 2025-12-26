/*
  File:             constants.js
  Description:      Shared constants used across multiple pages/modules.
  Author:           Adrian deCola
  Relative Path:    uSnapback/src/js/shared/constants.js
  Used by:          ../pages/*.js, other ../shared/*.js modules
*/

/* ----------------------------- */
/* ---------- Limits ----------- */
/* ----------------------------- */

// Max allowed nucleotides in the amplicon (excluding whitespace)
export const AMPLICON_LIMIT = 1000;

// Minimum amplicon length required by the app logic
export const MIN_AMP_LEN = 33;

// Minimum primer length (nt)
export const MIN_PRIMER_LEN = 12;

// Minimum number of bases BETWEEN forward and reverse primers
export const MIN_GAP_BETWEEN_PRIMERS = 7;

// SNV must be at least this many bases away from primer regions.
// (Your current rule is effectively a "≥ 3-bp gap" from primer binding sites.)
export const SNV_GAP = 3;

// Desired Tm bounds (°C)
export const TM_MIN = 40;
export const TM_MAX = 80;
