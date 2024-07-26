// Simple proxy server for testing IGV

// Proxy server from https://github.com/kasattejaswi/nodejs-proxy-server
// Copyright (c) 2021 Tejaswi Kasat

// Import of net module
const net = require("net")
const server = net.createServer()

server.on("connection", (clientToProxySocket) => {

    console.log("Client connected to proxy")

    clientToProxySocket.once("data", (data) => {

        console.log(`data -> ${data.toString()}`)

        let dataString = data.toString()
        let isTLSConnection = dataString.indexOf("CONNECT") !== -1

        let serverPort = 80
        let serverAddress
        if (isTLSConnection) {
            serverPort = 443
            serverAddress = dataString
                .split("CONNECT")[1]
                .split(" ")[1]
                .split(":")[0]
        } else {
            serverAddress = dataString.split("Host: ")[1].split("\r\n")[0]
        }
        console.log(serverAddress)

        // Require a password
        //if(dataString.indexOf("Proxy-Authorization") < 0) {
        //   clientToProxySocket.write("HTTP/1.1 407 Proxy requires authentication\r\n\r\n")
         //  return;
        //}

        // Creating a connection from proxy to destination server
        let proxyToServerSocket = net.createConnection(
            {
                host: serverAddress,
                port: serverPort,
            },
            () => {
                console.log("Proxy to server set up")
            }
        )


        if (isTLSConnection) {
            clientToProxySocket.write("HTTP/1.1 200 OK\r\n\r\n")
        } else {
            proxyToServerSocket.write(data)
        }

        clientToProxySocket.pipe(proxyToServerSocket)
        proxyToServerSocket.pipe(clientToProxySocket)


    })
})

server.on("error", (err) => {
    console.log("Some internal server error occurred")
    console.log(err)
})

server.on("close", () => {
    console.log("Client disconnected")
})

server.listen(9999)
