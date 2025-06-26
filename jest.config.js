// jest.config.js
import { createDefaultPreset } from 'ts-jest';

const tsJestTransformCfg = createDefaultPreset().transform;

/** @type {import("jest").Config} **/
export default {
	testEnvironment: 'jsdom', // We are testing browser code
	transform: {
		...tsJestTransformCfg,
	},
	setupFiles: ['./jest.setup.js'], // polyfill file
};
