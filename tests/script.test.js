// tests/script.test.js

import { JSDOM } from 'jsdom';

// Polyfill DOMParser globally, its used in parseTmFromResponse but its only
// availible globally in the browser and not in node (where Jest runs)
global.DOMParser = new JSDOM().window.DOMParser;

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

	// Constants
	SNV_BASE_BUFFER,
} from '../src/script.js';

/****************************************************************/
/*********************** Primary Function ***********************/
/****************************************************************/

/****************************************************************/
/********************* Secondary Functions **********************/
/****************************************************************/

describe('getStemTm()', () => {
	// Valid parameters
	test('returns 47.27 for sequence "gaaaaggagtgca" with no mismatch', async () => {
		const result = await getStemTm('GAAAAGGAGTGCA');
		expect(result).toBeCloseTo(47.27, 2);
	});

	test('returns 47.27 for sequence "gaaaaggagtgca" with null mismatch', async () => {
		const result = await getStemTm('GAAAAGGAGTGCA', null);
		expect(result).toBeCloseTo(47.27, 2);
	});

	test('returns 37.54 with valid mismatch', async () => {
		const mismatch = { position: 4, type: 'G' };
		const result = await getStemTm('GAAAAGGAGTGCA', mismatch);
		expect(result).toBeCloseTo(37.54, 2);
	});

	test('does NOT throw if mismatch is null', async () => {
		await expect(getStemTm('GAAAAGGAGTGCA', null)).resolves.not.toThrow();
	});

	// Invalid parameters
	test('throws if sequence is not a string', async () => {
		await expect(() => getStemTm(12345)).rejects.toThrow(
			'Invalid, empty, or non-string DNA sequence: "12345"'
		);
	});

	test('throws if sequence has invalid characters', async () => {
		await expect(() => getStemTm('GAXXXC')).rejects.toThrow(
			'Invalid, empty, or non-string DNA sequence: "GAXXXC"'
		);
	});

	test('throws if mismatch is not an object', async () => {
		await expect(() =>
			getStemTm('GAAAAGGAGTGC', 'not-an-object')
		).rejects.toThrow('Mismatch must be an object. Received: string');
	});

	test('throws if mismatch is an array', async () => {
		await expect(() =>
			getStemTm('GAAAAGGAGTGC', [{ position: 3, type: 'G' }])
		).rejects.toThrow('Mismatch must be an object. Received: object');
	});

	test('throws if mismatch is missing "position"', async () => {
		await expect(() =>
			getStemTm('GAAAAGGAGTGC', { type: 'G' })
		).rejects.toThrow(
			'Mismatch object missing required keys "position" and/or "type". Received: {"type":"G"}'
		);
	});

	test('throws if mismatch is missing "type"', async () => {
		await expect(() =>
			getStemTm('GAAAAGGAGTGC', { position: 3 })
		).rejects.toThrow(
			'Mismatch object missing required keys "position" and/or "type". Received: {"position":3}'
		);
	});

	test('throws if mismatch.position is out of bounds', async () => {
		await expect(() =>
			getStemTm('GAAAAGGAGTGC', { position: 100, type: 'G' })
		).rejects.toThrow(
			'Mismatch position 100 is invalid or out of bounds for sequence of length 12'
		);
	});

	test('throws if mismatch.position is negative', async () => {
		await expect(() =>
			getStemTm('GAAAAGGAGTGC', { position: -1, type: 'G' })
		).rejects.toThrow(
			'Mismatch position -1 is invalid or out of bounds for sequence of length 12'
		);
	});

	test('throws if mismatch.position is not an integer', async () => {
		await expect(() =>
			getStemTm('GAAAAGGAGTGC', { position: 2.5, type: 'G' })
		).rejects.toThrow(
			'Mismatch position 2.5 is invalid or out of bounds for sequence of length 12'
		);
	});

	test('throws if mismatch.type is invalid base', async () => {
		await expect(() =>
			getStemTm('GAAAAGGAGTGC', { position: 4, type: 'X' })
		).rejects.toThrow(
			'Mismatch type "X" must be one of "A", "T", "C", "G"'
		);
	});

	test('throws if mismatch.type is too long', async () => {
		await expect(() =>
			getStemTm('GAAAAGGAGTGC', { position: 4, type: 'GA' })
		).rejects.toThrow(
			'Mismatch type "GA" must be one of "A", "T", "C", "G"'
		);
	});

	test('throws if mismatch has extra keys', async () => {
		await expect(() =>
			getStemTm('GAAAAGGAGTGC', {
				position: 3,
				type: 'A',
				unexpected: true,
			})
		).rejects.toThrow(
			'Mismatch object contains unexpected key: "unexpected"'
		);
	});
});

/****************************************************************/
/*********************** Helper Functions ***********************/
/****************************************************************/

describe('snvTooCloseToPrimer()', () => {
	// Valid parameter behavior
	test('returns false if SNV is just far away enough from "left-side primer"', () => {
		expect(snvTooCloseToPrimer(5 + SNV_BASE_BUFFER, 5, 10, 50)).toBe(false);
		expect(snvTooCloseToPrimer(6 + SNV_BASE_BUFFER, 6, 10, 60)).toBe(false);
	});

	test('returns false when SNV is just far away enough from "right-side primer', () => {
		const seqLen = 50;
		const primerLen = 5;
		const compPrimerLen = 5;
		const snvIndex = seqLen - compPrimerLen - SNV_BASE_BUFFER - 1;
		expect(
			snvTooCloseToPrimer(snvIndex, primerLen, compPrimerLen, seqLen)
		).toBe(false);
	});

	test('returns false when SNV is comfortably away from both primers', () => {
		expect(snvTooCloseToPrimer(20, 5, 5, 50)).toBe(false);
		expect(snvTooCloseToPrimer(25, 10, 10, 70)).toBe(false);
	});

	test('returns true if SNV is too close to left primer', () => {
		expect(snvTooCloseToPrimer(7, 5, 10, 50)).toBe(true);
		expect(snvTooCloseToPrimer(4, 4, 5, 30)).toBe(true);
	});

	test('returns true if SNV is too close to right primer', () => {
		expect(snvTooCloseToPrimer(43, 5, 5, 50)).toBe(true);
	});

	test('returns true if SNV is both too close to left and right primer', () => {
		expect(snvTooCloseToPrimer(9, 9, 9, 20)).toBe(true);
	});

	test('returns true when SNV is exactly one base inside lower threshold', () => {
		const snvIndex = 5 + SNV_BASE_BUFFER - 1;
		expect(snvTooCloseToPrimer(snvIndex, 5, 5, 50)).toBe(true);
	});

	test('returns true when SNV is exactly one base outside upper threshold', () => {
		const snvIndex = 50 - 5 - SNV_BASE_BUFFER;
		expect(snvTooCloseToPrimer(snvIndex, 5, 5, 50)).toBe(true);
	});

	// Invalid input tests
	test('throws if any input is not a number', () => {
		expect(() => snvTooCloseToPrimer('1', 5, 5, 50)).toThrow(
			'snvIndex must be a finite number'
		);
		expect(() => snvTooCloseToPrimer(1, '5', 5, 50)).toThrow(
			'primerLen must be a finite number'
		);
		expect(() => snvTooCloseToPrimer(1, 5, null, 50)).toThrow(
			'compPrimerLen must be a finite number'
		);
		expect(() => snvTooCloseToPrimer(1, 5, 5, undefined)).toThrow(
			'seqLen must be a finite number'
		);
	});

	test('throws if any input is not an integer', () => {
		expect(() => snvTooCloseToPrimer(1.5, 5, 5, 50)).toThrow(
			'snvIndex must be an integer'
		);
		expect(() => snvTooCloseToPrimer(1, 5.2, 5, 50)).toThrow(
			'primerLen must be an integer'
		);
		expect(() => snvTooCloseToPrimer(1, 5, 5.9, 50)).toThrow(
			'compPrimerLen must be an integer'
		);
		expect(() => snvTooCloseToPrimer(1, 5, 5, 50.01)).toThrow(
			'seqLen must be an integer'
		);
	});

	test('throws for negative values', () => {
		expect(() => snvTooCloseToPrimer(-1, 5, 5, 50)).toThrow(
			'snvIndex must be non-negative'
		);
		expect(() => snvTooCloseToPrimer(1, -5, 5, 50)).toThrow(
			'primerLen must be non-negative'
		);
		expect(() => snvTooCloseToPrimer(1, 5, -5, 50)).toThrow(
			'compPrimerLen must be non-negative'
		);
		expect(() => snvTooCloseToPrimer(1, 5, 5, -50)).toThrow(
			'seqLen must be non-negative'
		);
	});

	test('throws if snvIndex is equal to seqLen', () => {
		expect(() => snvTooCloseToPrimer(50, 5, 5, 50)).toThrow(
			/snvIndex \(50\) cannot exceed or equal/
		);
	});

	test('throws if snvIndex > seqLen', () => {
		expect(() => snvTooCloseToPrimer(51, 5, 5, 50)).toThrow();
	});

	test('throws if total primer lengths exceed sequence length', () => {
		expect(() => snvTooCloseToPrimer(10, 25, 25, 50)).toThrow(
			/Primer lengths.*exceed or match/
		);
		expect(() => snvTooCloseToPrimer(10, 30, 30, 60)).toThrow();
	});

	// Edge cases (can change to throw errors if I want to implement differently)
	test('returns false when primers are 0-length and SNV is within bounds', () => {
		expect(snvTooCloseToPrimer(10, 0, 0, 20)).toBe(false);
		expect(snvTooCloseToPrimer(3, 0, 0, 10)).toBe(false);
	});

	test('returns true when primers are 0 and SNV is too close to sequence edge', () => {
		expect(snvTooCloseToPrimer(0, 0, 0, 10)).toBe(true);
		expect(snvTooCloseToPrimer(9, 0, 0, 10)).toBe(true);
	});
});

describe('buildMismatchSequenceForAPI()', () => {
	const shortSeq = 'GTCT'; // 4 bp
	const mediumSeq = 'CCTGAGGTTTTAACGGTTGTATGTAGGTGTATATTATCAG'; // 40 bp
	const longSeq =
		'TACTTTGGGATTACACCACCGTCTTTTCGTATGGAACAAGTCTGCCGCGAAACGTGCTGAGTGTACTAACTCGATCTTTCGAGCAAGGAGGTTTTCCGTCCCTGTCTTTTTGTAATCAAA'; // 120 bp
	const veryLongSeq =
		'GCTCGATATCACACCCTGGCCGCGTTCTTTAATTTCATATGAGGCGACTCGTTGCACCCCGTTAACAGAAATCTGAGTGTTATCCCGTTAAGCTAGCGCGTGTCTGCCCGGACACATTTCGCTCCCGAGCCACATTCCTCAAAATCCAACTCGGGGCGCGCGTTGTACTCCTAGACTCCAATATGTAAATACGTGCGTCA'; // 200 bp

	// ====== Valid cases ======
	test('correctly changes base for short sequence at position 0', () => {
		const mismatch = { position: 0, type: 'A' }; // insert T at index 0
		const result = buildMismatchSequenceForAPI(shortSeq, mismatch);
		expect(result.length).toBe(shortSeq.length);
		expect(result[0]).toBe('T');
		expect(result.slice(1)).toBe(shortSeq.slice(1));
	});

	test('correctly changes base for medium sequence in the middle', () => {
		const mismatch = { position: 20, type: 'T' }; // insert A at index 20
		const result = buildMismatchSequenceForAPI(mediumSeq, mismatch);
		expect(result.length).toBe(mediumSeq.length);
		expect(result[20]).toBe('A');
		expect(result.slice(0, 20)).toBe(mediumSeq.slice(0, 20));
		expect(result.slice(21)).toBe(mediumSeq.slice(21));
	});

	test('correctly changes base for long sequence near end', () => {
		const mismatch = { position: 119, type: 'C' }; // insert G at final index
		const result = buildMismatchSequenceForAPI(longSeq, mismatch);
		expect(result.length).toBe(longSeq.length);
		expect(result[119]).toBe('G');
		expect(result.slice(0, 119)).toBe(longSeq.slice(0, 119));
	});

	test('correctly changes base for very long sequence in middle', () => {
		const mismatch = { position: 100, type: 'G' }; // insert C at index 100
		const result = buildMismatchSequenceForAPI(veryLongSeq, mismatch);
		expect(result.length).toBe(veryLongSeq.length);
		expect(result[100]).toBe('C');
		expect(result.slice(0, 100)).toBe(veryLongSeq.slice(0, 100));
		expect(result.slice(101)).toBe(veryLongSeq.slice(101));
	});

	// ====== Invalid input cases ======
	test('throws on invalid base in sequence', () => {
		expect(() =>
			buildMismatchSequenceForAPI('ATBX', { position: 2, type: 'C' })
		).toThrow('Invalid input sequence');
	});

	test('throws for non-string sequence', () => {
		expect(() =>
			buildMismatchSequenceForAPI(null, { position: 0, type: 'A' })
		).toThrow('Invalid input sequence');
		expect(() =>
			buildMismatchSequenceForAPI(undefined, { position: 0, type: 'A' })
		).toThrow();
		expect(() =>
			buildMismatchSequenceForAPI(42, { position: 0, type: 'A' })
		).toThrow();
	});

	test('throws for missing mismatch keys', () => {
		expect(() => buildMismatchSequenceForAPI('ATCG', {})).toThrow(
			'Invalid mismatch object'
		);
		expect(() =>
			buildMismatchSequenceForAPI('ATCG', { position: 1 })
		).toThrow('Invalid mismatch object');
		expect(() =>
			buildMismatchSequenceForAPI('ATCG', { type: 'A' })
		).toThrow('Invalid mismatch object');
	});

	test('throws for out-of-bounds position', () => {
		expect(() =>
			buildMismatchSequenceForAPI('ATCG', { position: 4, type: 'A' })
		).toThrow('out of bounds');
		expect(() =>
			buildMismatchSequenceForAPI('ATCG', { position: -1, type: 'T' })
		).toThrow('out of bounds');
	});

	test('throws for invalid mismatch type values', () => {
		expect(() =>
			buildMismatchSequenceForAPI('ATCG', { position: 1, type: 'X' })
		).toThrow('Mismatch type "X" must be a single character');
		expect(() =>
			buildMismatchSequenceForAPI('ATCG', { position: 1, type: '' })
		).toThrow();
		expect(() =>
			buildMismatchSequenceForAPI('ATCG', { position: 1, type: 'AA' })
		).toThrow();
		expect(() =>
			buildMismatchSequenceForAPI('ATCG', { position: 1, type: 42 })
		).toThrow();
	});

	test('throws for lowercase mismatch type', () => {
		expect(() =>
			buildMismatchSequenceForAPI('ATCG', { position: 1, type: 'g' })
		).toThrow();
	});
});

describe('parseTmFromResponse()', () => {
	// Valid cases
	const validMatchHTML = `
		<html>
		<head></head>
		<body>
			<seq> gaaaaggagtgca </seq>
			<tm> 47.27 </tm>
		</body>
		</html>
	`;

	const validMismatchHTML = `
		<html>
		<head></head>
		<body>
			<seq> gaaaaggagtgca </seq>
			<tm> 47.04 </tm>
			<mmtm> 37.54 </mmtm>
			<warnings> Can only use Santalucia/Hicks parameter set with mismatches. Those parameters will be used here. </warnings>
		</body>
		</html>
	`;

	test('returns 47.27 when mismatch is not passed', () => {
		const result = parseTmFromResponse(validMatchHTML);
		expect(result).toBe(47.27);
	});

	test('returns 47.27 when mismatch is false', () => {
		const result = parseTmFromResponse(validMatchHTML, false);
		expect(result).toBe(47.27);
	});

	test('returns 37.54 when mismatch is true and <mmtm> is present', () => {
		const result = parseTmFromResponse(validMismatchHTML, true);
		expect(result).toBe(37.54);
	});

	// Invalid cases
	test('returns null when <tm> tag is missing and mismatch is false', () => {
		const html = `
			<html>
			<body>
				<seq> ATCG </seq>
			</body>
			</html>
		`;
		expect(parseTmFromResponse(html, false)).toBeNull();
	});

	test('returns null when <mmtm> tag is missing and mismatch is true', () => {
		const html = `
			<html>
			<body>
				<seq> ATCG </seq>
				<tm> 45.5 </tm>
			</body>
			</html>
		`;
		expect(parseTmFromResponse(html, true)).toBeNull();
	});

	test('returns null when <tm> is not a valid number', () => {
		const html = `
			<html><body><tm>not-a-number</tm></body></html>
		`;
		expect(parseTmFromResponse(html)).toBeNull();
	});

	test('returns null when <mmtm> is not a valid number', () => {
		const html = `
			<html><body><mmtm>NaN</mmtm></body></html>
		`;
		expect(parseTmFromResponse(html, true)).toBeNull();
	});

	test('returns null for empty string input', () => {
		expect(parseTmFromResponse('', true)).toBeNull();
		expect(parseTmFromResponse('', false)).toBeNull();
	});

	test('returns null for malformed HTML', () => {
		const html = `<<tm>>`;
		expect(parseTmFromResponse(html)).toBeNull();
	});
});

/****************************************************************/
/******************** DNA Utility Functions *********************/
/****************************************************************/

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
