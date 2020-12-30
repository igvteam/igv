

class IGVClient {

    constructor(host, port) {
        this.host = host || localhost;
        this.port = port || 60151;
    }

    async send(command, paramString) {

        var salt = Math.random(); // to prevent the browser from caching the response and preventing a relaunch if igv was shut down
        var localURL = `http://127.0.0.1:${this.port}/${command}?${decodeURIComponent(paramString)}&callback=igv.callBack();&salt=${salt}`;

    }

}