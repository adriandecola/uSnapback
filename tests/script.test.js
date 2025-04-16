/**
 * @jest-environment jsdom
 */

import {
	// Primary function
	createSnapback,

	// Secondary functions
	useTargetStrandsPrimerForComplement,
	evaluateSnapbackOptions,
	getStemTm,
	createStem,

	// Helper/logic functions
	snvTooCloseToPrimer,
	buildMismatchSequenceForAPI,
	parseTmFromResponse,

	// DNA utility functions
	isValidDNASequence,
	complementSequence,
	reverseComplement,
	revCompSNV,
} from '../src/script.js';

describe('isValidDNASequence()', () => {
	test('returns true for valid uppercase DNA sequences', () => {
		expect(isValidDNASequence('A')).toBe(true);
		expect(isValidDNASequence('ATCG')).toBe(true);
		expect(isValidDNASequence('GGGTTTAAAACCC')).toBe(true);
	});

	test('returns true for valid single uppercase characters', () => {
		expect(isValidDNASequence('A')).toBe(true);
		expect(isValidDNASequence('T')).toBe(true);
		expect(isValidDNASequence('C')).toBe(true);
		expect(isValidDNASequence('G')).toBe(true);
	});

	test('returns false for lowercase sequences', () => {
		expect(isValidDNASequence('a')).toBe(false);
		expect(isValidDNASequence('atcg')).toBe(false);
		expect(isValidDNASequence('gattaca')).toBe(false);
	});

	test('returns false for lowercase single characters', () => {
		expect(isValidDNASequence('a')).toBe(false);
		expect(isValidDNASequence('t')).toBe(false);
		expect(isValidDNASequence('c')).toBe(false);
		expect(isValidDNASequence('g')).toBe(false);
	});

	test('returns false for mixed-case sequences', () => {
		expect(isValidDNASequence('AtCg')).toBe(false);
		expect(isValidDNASequence('GgGgTtTt')).toBe(false);
	});

	test('returns false for sequences with numbers', () => {
		expect(isValidDNASequence('AT1G')).toBe(false);
		expect(isValidDNASequence('1234')).toBe(false);
		expect(isValidDNASequence('GATT4CA')).toBe(false);
	});

	test('returns false for sequences with special characters', () => {
		expect(isValidDNASequence('A*T^G')).toBe(false);
		expect(isValidDNASequence('@TCG')).toBe(false);
		expect(isValidDNASequence('ACGT!')).toBe(false);
		expect(isValidDNASequence('ACG$T')).toBe(false);
	});

	test('returns false for whitespace-containing strings', () => {
		expect(isValidDNASequence('AT CG')).toBe(false);
		expect(isValidDNASequence(' ATCG')).toBe(false);
		expect(isValidDNASequence('ATCG ')).toBe(false);
		expect(isValidDNASequence('AT\nCG')).toBe(false);
	});

	test('returns true for empty string (no invalid bases)', () => {
		expect(isValidDNASequence('')).toBe(true);
	});

	test('returns false for completely invalid strings', () => {
		expect(isValidDNASequence('xyz')).toBe(false);
		expect(isValidDNASequence('lmnop')).toBe(false);
		expect(isValidDNASequence('....')).toBe(false);
	});

	test('returns false for non-string input (if passed accidentally)', () => {
		expect(isValidDNASequence(1234)).toBe(false);
		expect(isValidDNASequence(null)).toBe(false);
		expect(isValidDNASequence(undefined)).toBe(false);
		expect(isValidDNASequence(true)).toBe(false);
		expect(isValidDNASequence(['A', 'T', 'C', 'G'])).toBe(false);
	});
});

describe('complementSequence()', () => {
	test('returns correct complement for simple valid sequences', () => {
		expect(complementSequence('A')).toBe('T');
		expect(complementSequence('T')).toBe('A');
		expect(complementSequence('C')).toBe('G');
		expect(complementSequence('G')).toBe('C');
		expect(complementSequence('ATCG')).toBe('TAGC');
		expect(complementSequence('GGGTTTAAAACCC')).toBe('CCCAAATTTTGGG');
	});

	test('throws error on lowercase input', () => {
		expect(() => complementSequence('atcg')).toThrow(
			'Invalid DNA sequence'
		);
		expect(() => complementSequence('a')).toThrow('Invalid DNA sequence');
	});

	test('throws error on mixed case', () => {
		expect(() => complementSequence('AtCg')).toThrow(
			'Invalid DNA sequence'
		);
	});

	test('throws error on invalid characters', () => {
		expect(() => complementSequence('AT1G')).toThrow(
			'Invalid DNA sequence'
		);
		expect(() => complementSequence('A*T^G')).toThrow(
			'Invalid DNA sequence'
		);
		expect(() => complementSequence('ACGT!')).toThrow(
			'Invalid DNA sequence'
		);
		expect(() => complementSequence('ACGT ')).toThrow(
			'Invalid DNA sequence'
		);
	});

	test('throws error on whitespace', () => {
		expect(() => complementSequence('A T G')).toThrow(
			'Invalid DNA sequence'
		);
		expect(() => complementSequence('AT\nCG')).toThrow(
			'Invalid DNA sequence'
		);
	});

	test('returns empty string if input is empty', () => {
		expect(complementSequence('')).toBe('');
	});

	test('throws error for non-string input', () => {
		expect(() => complementSequence(1234)).toThrow('Invalid DNA sequence');
		expect(() => complementSequence(null)).toThrow('Invalid DNA sequence');
		expect(() => complementSequence(undefined)).toThrow(
			'Invalid DNA sequence'
		);
		expect(() => complementSequence(['A', 'T'])).toThrow(
			'Invalid DNA sequence'
		);
	});
});
