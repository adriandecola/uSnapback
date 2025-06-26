// jest.setup.js

// A "polyfill" adds missing functionality to an environment — for example,
// browser APIs like TextEncoder, DOMParser, or fetch — so that Node.js can
// run browser-like code in tests. Typically imports are resolved before any
// code in a module/file runs (static imports)--the import statements are
// hoisted to the top of the fiel. Therefore, for modules that require other
// globals to be polyfilled, we must import them dynamically using the syntax
// seen below

// ------------------------------
// Polyfill TextEncoder/TextDecoder
// These are required by libraries like jsdom.
// Must be done before importing anything that relies on them.
// ------------------------------
import { TextEncoder, TextDecoder } from 'util';
global.TextEncoder = TextEncoder;
global.TextDecoder = TextDecoder;

// ------------------------------
// Polyfill DOMParser
// We defer this import so it only runs after TextEncoder is defined.
// This is done with the following dynamic import statement
// DOMParser is usually only available in browser environments.
// jsdom provides it here so we can simulate browser behavior.
// ------------------------------
const { JSDOM } = await import('jsdom');
global.DOMParser = new JSDOM().window.DOMParser;

// ------------------------------
// Polyfill fetch using node-fetch
// ------------------------------
import fetch from 'node-fetch';
if (!global.fetch) global.fetch = fetch;

// ------------------------------
// Set a global test timeout in Jest
// This prevents individual tests from running forever
// ------------------------------
import { jest } from '@jest/globals';
jest.setTimeout(10_000);
