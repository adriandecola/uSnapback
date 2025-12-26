/*
  File:             start.js
  Description:      Page logic for start.html
  Author:           Adrian deCola
  Relative Path:    uSnapback/src/js/pages/start.js
  Used by:          ../../pages/start.html
*/

import {
	AMPLICON_LIMIT,
	MIN_AMP_LEN,
	MIN_PRIMER_LEN,
	MIN_GAP_BETWEEN_PRIMERS,
	SNV_GAP,
	TM_MIN,
	TM_MAX,
} from '../shared/constants.js';

/*
  Note: using type="module" in the HTML script tag already avoids polluting global scope.
*/

const nextBtn = document.getElementById('nextBtn');

if (nextBtn) {
	nextBtn.addEventListener('click', () => {
		window.location.href = './amplicon.html';
	});
}
