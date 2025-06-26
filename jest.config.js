// jest.config.js
export default {
	testEnvironment: 'jsdom', // because I'm testing code that runs in the browser
	setupFiles: ['./jest.setup.js'], // injects globals that node 22 does not provide but jest expects
};
