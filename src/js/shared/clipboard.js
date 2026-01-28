/*
  File:             clipboard.js
  Description:      Clipboard helpers shared across pages
  Author:           Adrian deCola
  Relative Path:    uSnapback/src/js/shared/clipboard.js
*/

export function copyTextToClipboard(text) {
	if (!text) return Promise.resolve(false);

	// Prefer modern async clipboard API
	if (navigator.clipboard && navigator.clipboard.writeText) {
		return navigator.clipboard.writeText(text).then(
			() => true,
			(err) => {
				console.error('Clipboard copy failed:', err);
				return false;
			},
		);
	}

	// Fallback for older browsers
	const textarea = document.createElement('textarea');
	textarea.value = text;
	textarea.setAttribute('readonly', '');
	textarea.style.position = 'absolute';
	textarea.style.left = '-9999px';
	document.body.appendChild(textarea);
	textarea.select();

	let ok = false;
	try {
		ok = document.execCommand('copy');
	} catch (err) {
		console.error('execCommand copy failed:', err);
		ok = false;
	}

	document.body.removeChild(textarea);
	return Promise.resolve(ok);
}

export function wireCopyButton(buttonEl, targetEl, options = {}) {
	if (!buttonEl || !targetEl) return;
	const { statusEl, message = 'Copied to clipboard!', duration = 1500 } =
		options;

	buttonEl.addEventListener('click', () => {
		const text = (targetEl.textContent || '').trim();
		copyTextToClipboard(text).then((ok) => {
			if (!ok || !statusEl) return;

			statusEl.textContent = message;
			statusEl.classList.add('is-visible');

			if (statusEl._copyTimeout) {
				clearTimeout(statusEl._copyTimeout);
			}
			statusEl._copyTimeout = setTimeout(() => {
				statusEl.classList.remove('is-visible');
			}, duration);
		});
	});
}
