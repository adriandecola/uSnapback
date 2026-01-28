/*
  File:             clipboard.js
  Description:      Clipboard helpers shared across pages
  Author:           Adrian deCola
  Relative Path:    uSnapback/src/js/shared/clipboard.js
*/

export function copyTextToClipboard(text) {
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

export function wireCopyButton(buttonEl, targetEl) {
	if (!buttonEl || !targetEl) return;

	buttonEl.addEventListener('click', () => {
		const text = (targetEl.textContent || '').trim();
		copyTextToClipboard(text);
	});
}
