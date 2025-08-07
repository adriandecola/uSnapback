// build.js
import fs from 'fs';
import path from 'path';
import dotenv from 'dotenv';
import { fileURLToPath } from 'url';

// 1. Store enviroment variables to process.env
dotenv.config();

// 2. Read enviroment variaples
const API_URL = process.env.API_URL?.trim();
const PROXY_URL = process.env.PROXY_URL?.trim();
const USE_PROXY = (process.env.USE_PROXY || '').toLowerCase() === 'true';

// 3. Check correct enviroment variables exist
if (!API_URL) {
	console.error('Missing API_URL in .env');
	process.exit(1);
}
if (USE_PROXY && !PROXY_URL) {
	console.error('USE_PROXY=true but PROXY_URL is missing in .env');
	process.exit(1);
}

// 4. Get/create project-relative paths
const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const SRC_DIR = path.join(__dirname, 'src');
const DIST_DIR = path.join(__dirname, 'dist');

// 5. Clear existing dist/ folder if one exist
fs.rmSync(DIST_DIR, { recursive: true, force: true });
fs.mkdirSync(DIST_DIR, { recursive: true });

// 6. Copy all files in src/ to dist/, replacing enviroment variables
for (const file of fs.readdirSync(SRC_DIR)) {
	const srcPath = path.join(SRC_DIR, file);
	const distPath = path.join(DIST_DIR, file);

	if (fs.statSync(srcPath).isDirectory()) continue;

	let content = fs.readFileSync(srcPath, 'utf-8');

	// Replace __API_URL__ in script.js only
	if (file === 'script.js') {
		content = content
			.replace(/__API_URL__/g, API_URL)
			.replace(/__PROXY_URL__/g, PROXY_URL || '')
			.replace(/__USE_PROXY__/g, String(USE_PROXY)); // "true" | "false"
	}

	fs.writeFileSync(distPath, content);
}

console.log(`Build complete.`);
