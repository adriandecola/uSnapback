let stemResizeObserver = null;
let stemSnvInterval = null;

// Complementary base match and helper function
const COMP = { A: 'T', T: 'A', C: 'G', G: 'C' };

const complementBase = (b) =>
	b ? COMP[String(b).toUpperCase()] || String(b).toUpperCase() : '';
// Normalize base to uppercase string
const normBase = (b) => (b ? String(b).trim().toUpperCase() : '');
// Is it a correct DNA Base?
const isDnaBase = (b) => b === 'A' || b === 'C' || b === 'G' || b === 'T';

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
	stemResizeObserver =
		stemResizeObserver ||
		new ResizeObserver((entries) => {
			entries.forEach((entry) => updateStemLayout(entry.target));
		});
}

function clearStemSnvInterval() {
	if (stemSnvInterval) {
		clearInterval(stemSnvInterval);
		stemSnvInterval = null;
	}
}

stemSnvInterval = setInterval(() => {
	showWild = !showWild;
	applyState();
}, 2000);

export { renderStemDiagram };
