/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.oauth;

import biz.source_code.base64Coder.Base64Coder;
import org.broad.igv.Globals;
import org.broad.igv.batch.CommandExecutor;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.StringUtils;

import java.awt.*;
import java.io.*;
import java.net.ServerSocket;
import java.net.Socket;
import java.net.URLDecoder;
import java.nio.channels.ClosedByInterruptException;
import java.security.NoSuchAlgorithmException;
import java.util.*;

/**
 * NOT CURRENTLY USED -- developed for dynamic redirectURIs (dynamic random port numers), possible with Google but not Amazaon
 */

public class OauthListener implements Runnable {

    private static Logger log = LogManager.getLogger(OauthListener.class);

    private static final String CRLF = "\r\n";

    private int port;
    private OAuthProvider provider;
    private Thread listenerThread;


    public static synchronized void start(int port, OAuthProvider provider) {
        OauthListener listener = null;
        try {
            listener = new OauthListener(port, provider);
            listener.listenerThread.start();
        } catch (Exception e) {
            log.error(e);
        }
    }

    private OauthListener(int port, OAuthProvider provider) {
        this.port = port;
        this.provider = provider;
        listenerThread = new Thread(this);
    }

    public void run() {

        try (ServerSocket serverSocket = new ServerSocket(port);
             Socket clientSocket = serverSocket.accept();) {

            PrintWriter out = new PrintWriter(clientSocket.getOutputStream(), true);
            BufferedReader in = new BufferedReader(new InputStreamReader(clientSocket.getInputStream()));

            String inputLine = in.readLine();

            // Consume the remainder of the request, if any (typically request headers).   This is important to free the connection.
            String nextLine = in.readLine();
            while (nextLine != null && nextLine.length() > 0) {
                nextLine = in.readLine();
            }

            String[] tokens = inputLine.split(" ");
            if (tokens.length < 2) {
                sendTextResponse(out, "ERROR unexpected oauth request: " + inputLine);
            } else {

                String[] parts = tokens[1].split("\\?");
                Map<String, String> params = parts.length < 2 ? new HashMap() : parseParameters(parts[1]);

                if (params.containsKey("error")) {
                    sendTextResponse(out, "Error authorizing IGV: " + params.get("error"));
                } else if (params.containsKey("code")) {
                    provider.fetchAccessToken(params.get("code"));
                    sendTextResponse(out, "Authorization successful.  You may close this tab.");
                } else {
                    sendTextResponse(out, "Unsuccessful authorization response: " + inputLine);
                }
            }
        } catch (java.net.BindException e) {
            MessageUtils.showErrorMessage("Error opening listener for oAuth authorization", e);
        } catch (Exception e) {
            MessageUtils.showErrorMessage("Error opening listener for oAuth authorization", e);
            log.error(e);
        }
    }


    private static final String HTTP_RESPONSE = "HTTP/1.1 200 OK";
    private static final String HTTP_NO_RESPONSE = "HTTP/1.1 204 No Response";
    private static final String CONNECTION_CLOSE = "Connection: close";
    private static final String NO_CACHE = "Cache-Control: no-cache, no-store";
    private static final String ACCESS_CONTROL_ALLOW_ORIGIN = "Access-Control-Allow-Origin: *";

    private void sendTextResponse(PrintWriter out, String result) {
        sendHTTPResponse(out, result, "text/html", "GET");
    }

    private void sendHTTPResponse(PrintWriter out, String result, String contentType, String method) {

        out.print(result == null ? HTTP_NO_RESPONSE : HTTP_RESPONSE);
        out.print(CRLF);
        out.print(ACCESS_CONTROL_ALLOW_ORIGIN);
        out.print(CRLF);
        if (result != null) {
            out.print("Content-Type: " + contentType);
            out.print(CRLF);
            out.print("Content-Length: " + (result.length()));
            out.print(CRLF);
            out.print(NO_CACHE);
            out.print(CRLF);
            out.print(CONNECTION_CLOSE);
            out.print(CRLF);

            if (!method.equals("HEAD")) {
                out.print(CRLF);
                out.print(result);
                out.print(CRLF);
            }
        } else {
            out.print(CRLF);
        }
        out.close();
    }


    /**
     * Parse the html parameter string into a set of key-value pairs.  Parameter values are
     * url decoded with the exception of the "locus" parameter.
     *
     * @param parameterString
     * @return
     */
    private Map<String, String> parseParameters(String parameterString) {

        // Do a partial decoding now (ampersands only)
        parameterString = parameterString.replace("&amp;", "&");

        HashMap<String, String> params = new HashMap();
        String[] kvPairs = parameterString.split("&");
        for (String kvString : kvPairs) {
            // Split on the first "=",  all others are part of the parameter value
            String[] kv = kvString.split("=", 2);
            if (kv.length == 1) {
                params.put(kv[0], null);
            } else {
                String key = StringUtils.decodeURL(kv[0]);
                String value = StringUtils.decodeURL(kv[1]);
                params.put(key, value);
            }
        }
        return params;
    }

    /**
     * Find and available port*
     * @return
     */
    public static int findFreePort() {
        for (int port = 49152; port < 65536; port++) {
            try (ServerSocket serverSocket = new ServerSocket(port)) {
                if (serverSocket != null && serverSocket.getLocalPort() == port) {
                    return port;
                }
            } catch (IOException  e) {
                //Port is not available, try again
            }
        }
        return -1;
    }

}
