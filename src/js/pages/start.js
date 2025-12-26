/*
  File:             start.js
  Description:      Page logic for start.html
  Author:           Adrian deCola
  Relative Path:    uSnapback/src/js/pages/start.js
  Used by:          ../../pages/start.html
*/

// Don't need to import any constants for this simple page

/*
  Note: using type="module" in the HTML script tag already avoids polluting global scope.
*/

const nextBtn = document.getElementById('nextBtn');

if (nextBtn) {
	nextBtn.addEventListener('click', () => {
		window.location.href = './amplicon.html';
	});
}
