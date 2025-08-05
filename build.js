// build.js
import fs from 'fs';
import path from 'path';
import dotenv from 'dotenv';
import { fileURLToPath } from 'url';

dotenv.config();

const API_URL = process.env.API_URL;
if (!API_URL) {
	console.error('Missing API_URL in .env');
	process.exit(1);
}

// Get/create project-relative paths
const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const SRC_DIR = path.join(__dirname, 'src');
const DIST_DIR = path.join(__dirname, 'dist');

// Clear existing dist/ folder if one exist
fs.rmSync(DIST_DIR, { recursive: true, force: true });
fs.mkdirSync(DIST_DIR, { recursive: true });

// Copy all files in src/ to dist/, replacing API endpoints
for (const file of fs.readdirSync(SRC_DIR)) {
	const srcPath = path.join(SRC_DIR, file);
	const distPath = path.join(DIST_DIR, file);

	if (fs.statSync(srcPath).isDirectory()) continue;

	let content = fs.readFileSync(srcPath, 'utf-8');

	// Replace __API_URL__ in script.js only
	if (file === 'script.js') {
		content = content.replace(/__API_URL__/g, API_URL);
	}

	fs.writeFileSync(distPath, content);
}

console.log(`🎉 Build complete. Injected API_URL: ${API_URL}`);
