// vite.config.js
import { defineConfig } from 'vite';

export default defineConfig(() => {
	const isTesting = process.env.TEST === 'true';

	return {
		root: 'dist',
		server: {
			port: 8000,
			open: isTesting ? false : '/pages/start.html',
		},
	};
});
