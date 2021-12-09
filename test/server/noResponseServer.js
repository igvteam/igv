/**
 * Simulates a server that does not return a response.
 *
 * To run
 *     node noResponseServer.js
 *
 * URL :  http://localhost:9999
 *
 * @type {module:http}
 */

const http = require("http")

const host = 'localhost'
const port = 9999

const server = http.createServer((req, res) => {
    //res.writeHead(200)
   // res.end("Response from server")
})
server.listen(port, host, () => {
    console.log(`Server is running on http://${host}:${port}`)
})

