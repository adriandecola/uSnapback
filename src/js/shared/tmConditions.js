/*
  File:             tmConditions.js
  Description:      Helpers for reading saved Tm calculation conditions
*/

import {
	DEFAULT_MAGNESIUM_MM,
	DEFAULT_MONOVALENT_MM,
} from './constants.js';

function readStoredConcentration(storage, key, fallback) {
	const raw = storage.getItem(key);
	if (raw == null || raw.trim() === '') return fallback;

	const value = Number(raw);
	return Number.isFinite(value) && value >= 0 ? value : fallback;
}

/**
 * Reads the last valid ion settings. Before the user visits the Tm page,
 * the requested defaults are used for primer previews.
 */
export function readStoredTmConditions(storage) {
	if (!storage || typeof storage.getItem !== 'function') {
		throw new Error('A storage object with getItem() is required.');
	}

	return {
		magnesiumMm: readStoredConcentration(
			storage,
			'magnesiumMm',
			DEFAULT_MAGNESIUM_MM,
		),
		monovalentMm: readStoredConcentration(
			storage,
			'monovalentMm',
			DEFAULT_MONOVALENT_MM,
		),
	};
}
