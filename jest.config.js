// jest.config.js
/** @type {import('ts-jest').JestConfigWithTsJest} */
export default {
	preset: 'ts-jest/presets/default-esm',
	testEnvironment: 'jsdom',
	extensionsToTreatAsEsm: ['.js'],
	globals: {
		'ts-jest': {
			useESM: true,
		},
	},
	moduleNameMapper: {
		'^(\\.{1,2}/.*)\\.js$': '$1',
	},
};
