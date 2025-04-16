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

describe('reverseComplement()', () => {
	// Simple valid tests
	test('returns correct reverse complement for valid sequences', () => {
		expect(reverseComplement('A')).toBe('T');
		expect(reverseComplement('T')).toBe('A');
		expect(reverseComplement('C')).toBe('G');
		expect(reverseComplement('G')).toBe('C');
		expect(reverseComplement('ATCG')).toBe('CGAT');
		expect(reverseComplement('GGGTTTAAA')).toBe('TTTAAACCC');
		expect(reverseComplement('ATGCG')).toBe('CGCAT');
		expect(reverseComplement('TACGT')).toBe('ACGTA');
		expect(reverseComplement('CCGGA')).toBe('TCCGG');
		expect(reverseComplement('TTTGG')).toBe('CCAAA');
	});

	// Long sequences
	test('handles long valid sequences', () => {
		const input = 'ATCGGGTTAACCGGATTCGAGCTTAAACCGGGTTTAAACC';
		const expected = 'GGTTTAAACCCGGTTTAAGCTCGAATCCGGTTAACCCGAT';
		expect(reverseComplement(input)).toBe(expected);
	});

	// Edge case: empty input
	test('returns empty string for empty input', () => {
		expect(reverseComplement('')).toBe('');
	});

	// Case sensitivity (should throw on non-uppercase)
	test('throws error on lowercase or mixed-case inputs', () => {
		expect(() => reverseComplement('atcg')).toThrow('Invalid DNA sequence');
		expect(() => reverseComplement('a')).toThrow('Invalid DNA sequence');
		expect(() => reverseComplement('AtCg')).toThrow('Invalid DNA sequence');
		expect(() => reverseComplement('GgGg')).toThrow('Invalid DNA sequence');
	});

	// Invalid characters
	test('throws error on invalid characters or formatting', () => {
		expect(() => reverseComplement('ATXG')).toThrow('Invalid DNA sequence');
		expect(() => reverseComplement('A T')).toThrow('Invalid DNA sequence');
		expect(() => reverseComplement('ACG$')).toThrow('Invalid DNA sequence');
		expect(() => reverseComplement('1234')).toThrow('Invalid DNA sequence');
	});

	// Non-string input
	test('throws error for non-string input', () => {
		expect(() => reverseComplement(42)).toThrow('Invalid DNA sequence');
		expect(() => reverseComplement(null)).toThrow('Invalid DNA sequence');
		expect(() => reverseComplement(undefined)).toThrow(
			'Invalid DNA sequence'
		);
		expect(() => reverseComplement(['A', 'T', 'C', 'G'])).toThrow(
			'Invalid DNA sequence'
		);
	});

	// Palindrome
	test('reverseComplement(reverseComplement(seq)) === seq', () => {
		const sequences = [
			'A',
			'ATCG',
			'GGGTTTAAA',
			'CCGGAATTC',
			'ATCGGGTTAACCGGATTCGAGCTTAAACCGGGTTTAAACC',
		];
		for (const seq of sequences) {
			expect(reverseComplement(reverseComplement(seq))).toBe(seq);
		}
	});
});

describe('revCompSNV()', () => {
	// Valid SNV inputs — variety of lengths
	test('returns correct reverse complement for SNVs in various length sequences', () => {
		expect(revCompSNV({ index: 0, variantBase: 'A' }, 5)).toEqual({
			index: 4,
			variantBase: 'T',
		});
		expect(revCompSNV({ index: 1, variantBase: 'T' }, 6)).toEqual({
			index: 4,
			variantBase: 'A',
		});
		expect(revCompSNV({ index: 2, variantBase: 'C' }, 7)).toEqual({
			index: 4,
			variantBase: 'G',
		});
		expect(revCompSNV({ index: 3, variantBase: 'G' }, 8)).toEqual({
			index: 4,
			variantBase: 'C',
		});
		expect(revCompSNV({ index: 4, variantBase: 'A' }, 9)).toEqual({
			index: 4,
			variantBase: 'T',
		});
		expect(revCompSNV({ index: 5, variantBase: 'T' }, 11)).toEqual({
			index: 5,
			variantBase: 'A',
		});
		expect(revCompSNV({ index: 6, variantBase: 'C' }, 12)).toEqual({
			index: 5,
			variantBase: 'G',
		});
		expect(revCompSNV({ index: 7, variantBase: 'G' }, 13)).toEqual({
			index: 5,
			variantBase: 'C',
		});
		expect(revCompSNV({ index: 8, variantBase: 'A' }, 14)).toEqual({
			index: 5,
			variantBase: 'T',
		});
		expect(revCompSNV({ index: 9, variantBase: 'T' }, 15)).toEqual({
			index: 5,
			variantBase: 'A',
		});
	});

	// Valid SNV inputs
	test('returns correct reverse complement for valid SNVs in 10-length sequence', () => {
		expect(revCompSNV({ index: 0, variantBase: 'A' }, 10)).toEqual({
			index: 9,
			variantBase: 'T',
		});
		expect(revCompSNV({ index: 1, variantBase: 'A' }, 10)).toEqual({
			index: 8,
			variantBase: 'T',
		});
		expect(revCompSNV({ index: 5, variantBase: 'A' }, 10)).toEqual({
			index: 4,
			variantBase: 'T',
		});
		expect(revCompSNV({ index: 9, variantBase: 'A' }, 10)).toEqual({
			index: 0,
			variantBase: 'T',
		});
		expect(revCompSNV({ index: 0, variantBase: 'T' }, 10)).toEqual({
			index: 9,
			variantBase: 'A',
		});
		expect(revCompSNV({ index: 5, variantBase: 'C' }, 10)).toEqual({
			index: 4,
			variantBase: 'G',
		});
		expect(revCompSNV({ index: 7, variantBase: 'G' }, 10)).toEqual({
			index: 2,
			variantBase: 'C',
		});
	});

	// Edge case: index in the middle of even and odd length seqs
	test('returns correct for center index', () => {
		expect(revCompSNV({ index: 4, variantBase: 'G' }, 9)).toEqual({
			index: 4,
			variantBase: 'C',
		});
		expect(revCompSNV({ index: 5, variantBase: 'C' }, 10)).toEqual({
			index: 4,
			variantBase: 'G',
		});
	});

	// Very short sequences
	test('works on shortest valid sequences (length 1-4)', () => {
		expect(revCompSNV({ index: 0, variantBase: 'A' }, 1)).toEqual({
			index: 0,
			variantBase: 'T',
		});
		expect(revCompSNV({ index: 1, variantBase: 'T' }, 3)).toEqual({
			index: 1,
			variantBase: 'A',
		});
		expect(revCompSNV({ index: 2, variantBase: 'G' }, 4)).toEqual({
			index: 1,
			variantBase: 'C',
		});
	});

	// Long sequence case
	test('works on long sequence', () => {
		expect(revCompSNV({ index: 49, variantBase: 'A' }, 100)).toEqual({
			index: 50,
			variantBase: 'T',
		});
		expect(revCompSNV({ index: 0, variantBase: 'G' }, 100)).toEqual({
			index: 99,
			variantBase: 'C',
		});
		expect(revCompSNV({ index: 99, variantBase: 'C' }, 100)).toEqual({
			index: 0,
			variantBase: 'G',
		});
	});

	// Invalid SNV objects
	test('throws for non-object snvSite input', () => {
		expect(() => revCompSNV(null, 10)).toThrow('Invalid SNV object');
		expect(() => revCompSNV(undefined, 10)).toThrow('Invalid SNV object');
		expect(() => revCompSNV([], 10)).toThrow('Invalid SNV object');
		expect(() => revCompSNV(42, 10)).toThrow('Invalid SNV object');
	});

	test('throws if snvSite is missing keys or has extra keys', () => {
		expect(() => revCompSNV({ index: 5 }, 10)).toThrow(
			'Invalid SNV object'
		);
		expect(() => revCompSNV({ variantBase: 'A' }, 10)).toThrow(
			'Invalid SNV object'
		);
		expect(() =>
			revCompSNV({ index: 5, variantBase: 'A', extra: true }, 10)
		).toThrow('Invalid SNV object');
	});

	// Invalid index
	test('throws for invalid index values', () => {
		expect(() => revCompSNV({ index: -1, variantBase: 'A' }, 10)).toThrow(
			'Invalid SNV index'
		);
		expect(() => revCompSNV({ index: 10, variantBase: 'T' }, 10)).toThrow(
			'Invalid SNV index'
		);
		expect(() => revCompSNV({ index: 5.5, variantBase: 'C' }, 10)).toThrow(
			'Invalid SNV index'
		);
		expect(() => revCompSNV({ index: '5', variantBase: 'G' }, 10)).toThrow(
			'Invalid SNV index'
		);
	});

	// Invalid base
	test('throws for invalid variantBase', () => {
		expect(() => revCompSNV({ index: 5, variantBase: 'X' }, 10)).toThrow(
			'Invalid variant base'
		);
		expect(() => revCompSNV({ index: 5, variantBase: 'g' }, 10)).toThrow(
			'Invalid variant base'
		);
		expect(() => revCompSNV({ index: 5, variantBase: '' }, 10)).toThrow(
			'Invalid variant base'
		);
		expect(() => revCompSNV({ index: 5, variantBase: 'AA' }, 10)).toThrow(
			'Invalid variant base'
		);
		expect(() => revCompSNV({ index: 5, variantBase: 1 }, 10)).toThrow(
			'Invalid variant base'
		);
	});
});
