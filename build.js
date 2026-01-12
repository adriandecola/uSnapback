// build.js
import fs from 'fs';
import path from 'path';
import dotenv from 'dotenv';
import { fileURLToPath } from 'url';

// 1. Store enviroment variables to process.env
dotenv.config();

// 2. Read enviroment variables
const DEPLOY_USNAPBACK =
	(process.env.DEPLOY_USNAPBACK || '').toLowerCase() === 'true';

const API_URL_UTAH =
	process.env.API_URL_UTAH?.trim() || process.env.API_URL?.trim();
const API_URL_USNAPBACK = process.env.API_URL_USNAPBACK?.trim();

const API_URL = DEPLOY_USNAPBACK ? API_URL_USNAPBACK : API_URL_UTAH;

const PROXY_URL = process.env.PROXY_URL?.trim();
const USE_PROXY = (process.env.USE_PROXY || '').toLowerCase() === 'true';

const USE_TOKEN = (process.env.USE_TOKEN || '').toLowerCase() === 'true';
const API_TOKEN = process.env.API_TOKEN?.trim() || '';

// 3. Check correct enviroment variables exist
if (!API_URL) {
	console.error(
		DEPLOY_USNAPBACK
			? 'DEPLOY_USNAPBACK=true but API_URL_USNAPBACK is missing in .env'
			: 'DEPLOY_USNAPBACK=false but API_URL_UTAH (or API_URL) is missing in .env'
	);
	process.exit(1);
}

if (USE_TOKEN && !API_TOKEN) {
	console.error('USE_TOKEN=true but API_TOKEN is missing in .env');
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

// 6. Copy everything from src/ into dist/ (including folders).
//    Only exception: in src/script.js we replace __API_URL__ / __PROXY_URL__ / __USE_PROXY__.

const foldersToVisit = ['']; // we start at /src

while (foldersToVisit.length > 0) {
	const relFolder = foldersToVisit.pop(); // e.g. '', 'pages', 'styles'
	const srcFolder = path.join(SRC_DIR, relFolder); // absolute path inside SRC_DIR

	// Read this folder and distinguish files vs subfolders
	const entries = fs.readdirSync(srcFolder, { withFileTypes: true });

	for (const entry of entries) {
		// macOS sometimes adds this file everywhere; ignore it
		if (entry.name === '.DS_Store') continue;

		// Build path for this entry, relative to src/
		const relPath = path.join(relFolder, entry.name); // e.g. 'pages/start.html'
		const srcPath = path.join(SRC_DIR, relPath);
		const distPath = path.join(DIST_DIR, relPath);

		// If it’s a folder:
		//   1) create the matching folder in dist/
		//   2) remember to walk into it
		if (entry.isDirectory()) {
			fs.mkdirSync(distPath, { recursive: true });
			foldersToVisit.push(relPath);
			continue;
		}

		// If it’s a file, make sure dist’s parent folder exists before writing/copying
		fs.mkdirSync(path.dirname(distPath), { recursive: true });

		// Only rewrite placeholders in the top-level script.js (src/script.js)
		// Everything else is copied byte-for-byte.
		const normalizedRelPath = relPath.replaceAll('\\', '/');
		if (normalizedRelPath === 'script.js') {
			let content = fs.readFileSync(srcPath, 'utf-8');

			content = content
				.replace(/__API_URL__/g, JSON.stringify(API_URL))
				.replace(/__PROXY_URL__/g, JSON.stringify(PROXY_URL || ''))
				.replace(/__USE_PROXY__/g, USE_PROXY ? 'true' : 'false')
				.replace(/__USE_TOKEN__/g, USE_TOKEN ? 'true' : 'false')
				.replace(/__API_TOKEN__/g, JSON.stringify(API_TOKEN));

			fs.writeFileSync(distPath, content, 'utf-8');
		} else {
			// Copy the file exactly as-is
			fs.copyFileSync(srcPath, distPath);
		}
	}
}

// 7. Complete
console.log(`Build complete.`);
