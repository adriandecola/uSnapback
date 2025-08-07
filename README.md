# uSnapback
An application to design a snapback primer for single base variant asymmetrical PCR. 

To test script.js, comment out the CORs proxy script in script.js and then run `npm run test` or `npm test`.
* This part is currently not working as I can't get those requests to come from port 8000 *

To build the dist/ folder and corresponding distribution files:
Ensure that you have a .env file with the proper key API_URL. This should be the API endpoint that you will hit, whether a middle man or not. Make sure this file is in the top level directory (same as src/ and dist/)
Then run `npm build`. 

To test the frontend/browser application make sure you are on port 8000. (Have not yet tested this)
