# uSnapback

An application to design a snapback primer for single base variant asymmetrical PCR.

To test script.js, via script.test.js, using Jest, comment out the CORs proxy script in script.js and then run `npm run test` or `npm test`. This builds the files, serves the application on port 8000, and then runs the test suite.

To build the dist/ folder and corresponding distribution files:

-   Ensure that you have a .env file with the proper key API_URL.
-   If you need to use a proxy make sure you also have USE_PROXY and PROXY_URL in the .env file. If no proxy is used, API_URL will be the endpoint you that is hit with the proper query. If you do use a proxy the URL hit will be PROXY_URL with the query 'url' with the value of the encoded API_URL target/query.
-   Make sure the .env file is in the top level directory (same as src/ and dist/)
    Then run `npm build`.

To test the frontend/browser application you can run `npm run serve` which serves the application on port 8000. Remember that this DOES NOT REBUILD your application automatically and will serve the start.html file in the dist/ directory.
