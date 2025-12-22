/*
  File:             start.js
  Description:      Page logic for start.html
  Author:           Adrian deCola
  Relative Path:    uSnapback/src/scripts/start.js
  Linked files:     ../pages/start.html
*/

/*
  Note: using type="module" in the HTML already avoids polluting global scope.
  I still keep the IIFE style you’ve been using.
*/

const nextBtn = document.getElementById('nextBtn');

if (nextBtn) {
	nextBtn.addEventListener('click', () => {
		window.location.href = './amplicon.html';
	});
}
