/*
  File:             constants.js
  Description:      Shared constants for uSnapback pages
  Author:           Adrian deCola
  Relative Path:    uSnapback/src/js/shared/constants.js
*/

export const AMPLICON_LIMIT = 1000;

/* Minimum amplicon length accepted by the app */
export const MIN_AMP_LEN = 33;

/* Minimum primer length (nt) */
export const MIN_PRIMER_LEN = 12;

/* Minimum number of bases BETWEEN forward and reverse primer sites */
export const MIN_GAP_BETWEEN_PRIMERS = 7;

/*
  Minimum number of bases between the SNV and a primer-binding site.
  This matches your current logic:
  - idx < fwdLen + 3  => too close to forward primer
  - idx > seqLen - revLen - 4 => too close to reverse primer
*/
export const SNV_GAP = 3;

/* Desired Tm range (°C) */
export const TM_MIN = 40;
export const TM_MAX = 80;

/* Default ionic conditions for snapback Tm estimates (mM): free Mg²⁺ and
   total monovalent-cation concentration, respectively. */
export const DEFAULT_MAGNESIUM_MM = 3.0;
export const DEFAULT_MONOVALENT_MM = 13.7;

/* Allowed DNA bases */
export const BASES = ['A', 'C', 'G', 'T'];
