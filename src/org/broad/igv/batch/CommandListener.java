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

package org.broad.igv.batch;

import biz.source_code.base64Coder.Base64Coder;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.ga4gh.OAuthUtils;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.StringUtils;

import java.awt.*;
import java.io.*;
import java.net.ServerSocket;
import java.net.Socket;
import java.net.URLDecoder;
import java.nio.channels.ClosedByInterruptException;
import java.security.NoSuchAlgorithmException;
import java.util.*;

public class CommandListener implements Runnable {

    public static final String OK = "OK";

    public static int currentListenerPort = -1;

    private static Logger log = Logger.getLogger(CommandListener.class);

    private static CommandListener listener;
    private static final String CRLF = "\r\n";

    private int port = -1;
    private ServerSocket serverSocket = null;
    private Socket clientSocket = null;
    private Thread listenerThread;
    boolean halt = false;


    /**
     * Different keys which can be used to specify a file to load
     */
    public static Set<String> fileParams;
    public static Set<String> indexParams;

    static {
        String[] fps = new String[]{"file", "bigDataURL", "sessionURL", "dataURL"};
        fileParams = new LinkedHashSet<String>(Arrays.asList(fps));
        fileParams = Collections.unmodifiableSet(fileParams);

        indexParams = new HashSet<String>(Arrays.asList("index"));
    }

    public static synchronized void start(int port) {
        listener = new CommandListener(port);
        listener.listenerThread.start();
    }


    public static synchronized void halt() {
        if (listener != null) {
            listener.halt = true;
            listener.listenerThread.interrupt();
            listener.closeSockets();
            listener = null;
        }
    }

    private CommandListener(int port) {
        this.port = port;
        listenerThread = new Thread(this);
    }

    /**
     * Loop forever, processing client requests synchronously.  The server is single threaded.
     */
    public void run() {

        CommandExecutor cmdExe = new CommandExecutor();

        try {
            serverSocket = new ServerSocket(port);
            log.info("Listening on port " + port);
            currentListenerPort = port;

            while (!halt) {
                clientSocket = serverSocket.accept();
                processClientSession(cmdExe);
                if (clientSocket != null) {
                    try {
                        clientSocket.close();
                        clientSocket = null;
                    } catch (IOException e) {
                        log.error("Error in client socket loop", e);
                    }
                }
            }


        } catch (java.net.BindException e) {
            log.error(e);
            currentListenerPort = -1;
        } catch (ClosedByInterruptException e) {
            log.error(e);

        } catch (IOException e) {
            if (!halt) {
                log.error("IO Error on port socket ", e);
            }
        }
    }

    /**
     * Process a client session.  Loop continuously until client sends the "halt" message, or closes the connection.
     *
     * @param cmdExe
     * @throws IOException
     */
    private void processClientSession(CommandExecutor cmdExe) throws IOException {
        PrintWriter out = null;
        BufferedReader in = null;
        try {
            out = new PrintWriter(clientSocket.getOutputStream(), true);
            in = new BufferedReader(new InputStreamReader(clientSocket.getInputStream()));
            String inputLine;


            while (!halt && (inputLine = in.readLine()) != null) {

                String cmd = inputLine;
                if (cmd.startsWith("GET")) {

                    // Consume the remainder of the request, if any.   This is important to free the connection.
                    Map<String, String> headers = new HashMap<String, String>();
                    String nextLine = in.readLine();
                    while (nextLine != null && nextLine.length() > 0) {
                        nextLine = in.readLine();
                        String[] tokens = Globals.colonPattern.split(nextLine, 2);
                        if (tokens.length == 2) {
                            headers.put(tokens[0].trim(), tokens[1].trim());
                        }
                    }

                    log.info(cmd);

                    String command = null;
                    Map<String, String> params = null;
                    String[] tokens = inputLine.split(" ");
                    if (tokens.length < 2) {
                        sendTextResponse(out, "ERROR unexpected command line: " + inputLine);
                        return;
                    } else {
                        String[] parts = tokens[1].split("\\?");
                        command = parts[0];
                        params = parts.length < 2 ? new HashMap() : parseParameters(parts[1]);
                    }

                    // Detect google oauth callback
                    if (command.equals("/oauthCallback")) {
                        if (params.containsKey("code")) {
                            OAuthUtils.getInstance().setAuthorizationCode(params.get("code"));
                        } else if (params.containsKey("token")) {
                            OAuthUtils.getInstance().setAccessToken(params.get("token"));
                        }
                        sendTextResponse(out, "OK");
                    } else {

                        // If a callback (javascript) function is specified write it back immediately.  This function
                        // is used to cancel a timeout handler
                        String callback = params.get("callback");
                        if (callback != null) {
                            sendJavascriptResponse(out, callback);
                        }

                        // Process the request.


                        String result = processGet(command, params, cmdExe);

                        // If no callback was specified write back response now
                        if (callback == null) {
                            // We send no response if result is "ok".
                            if (result.equals(OK)) result = null;
                            sendTextResponse(out, result);
                        }
                    }

                    // http sockets are used for one request only
                    return;
                    // }

                } else {
                    // Port command
                    Globals.setBatch(true);
                    Globals.setSuppressMessages(true);
                    final String response = cmdExe.execute(inputLine);
                    out.println(response);
                    out.flush();
                }
            }
        } catch (IOException e) {
            log.error("Error processing client session", e);
        } finally {
            Globals.setSuppressMessages(false);
            Globals.setBatch(false);
            if (out != null) out.close();
            if (in != null) in.close();
        }
    }


    private void closeSockets() {
        if (clientSocket != null) {
            try {
                clientSocket.close();
                clientSocket = null;
            } catch (IOException e) {
                log.error("Error closing clientSocket", e);
            }
        }

        if (serverSocket != null) {
            try {
                serverSocket.close();
                serverSocket = null;
            } catch (IOException e) {
                log.error("Error closing server socket", e);
            }
        }
    }


    private static final String HTTP_RESPONSE = "HTTP/1.1 200 OK";
    private static final String HTTP_NO_RESPONSE = "HTTP/1.1 204 No Response";
    private static final String CONNECTION_CLOSE = "Connection: close";
    private static final String NO_CACHE = "Cache-Control: no-cache, no-store";
    private static final String ACCESS_CONTROL_ALLOW_ORIGIN = "Access-Control-Allow-Origin: *";

    private void sendJavascriptResponse(PrintWriter out, String result) {
        sendHTTPResponse(out, result, "application/javascript");

    }

    private void sendTextResponse(PrintWriter out, String result) {
        sendHTTPResponse(out, result, "text/html");
    }

    private void sendHTTPResponse(PrintWriter out, String result, String contentType) {

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
            out.print(CRLF);
            out.print(result);
            out.print(CRLF);
        }
        out.close();
    }

    /**
     * Process an http get request.
     */

    private String processGet(String command, Map<String, String> params, CommandExecutor cmdExe) throws IOException {

        String result = OK;
        final Frame mainFrame = IGV.getMainFrame();

        // Trick to force window to front, the setAlwaysOnTop works on a Mac,  toFront() does nothing.
        mainFrame.toFront();
        mainFrame.setAlwaysOnTop(true);
        mainFrame.setAlwaysOnTop(false);

        //track%20type=bigBed%20name=%27hES_HUES1_p28.RRBS_CpG_meth%27%20description=%27RRBS%20CpG%20methylation%20for%20hES_HUES1_p28.RRBS%27%20visibility=4%20useScore=1%20color=0,60,120

        /** from what server was IGV started? Used to link to other tools/apps */

        String server = params.get("server");
        if (server == null || server.trim().length() < 1)
            server = PreferenceManager.getInstance().get(PreferenceManager.IONTORRENT_SERVER);
        else {
            PreferenceManager.getInstance().put(PreferenceManager.IONTORRENT_SERVER, server);
        }
        if (command.equals("/load")) {
            String file = null;
            for (String fp : fileParams) {
                file = params.get(fp);
                if (file != null) break;
            }

            String genome = params.get("genome");
            if (genome == null) {
                genome = params.get("db");  // <- UCSC track line param
            }

            if (genome != null) {
                IGV.getInstance().loadGenomeById(genome);
            }

            if (file != null) {
                PreferenceManager.getInstance().put(PreferenceManager.IONTORRENT_RESULTS, file);


                String mergeValue = params.get("merge");
                if (mergeValue != null) mergeValue = URLDecoder.decode(mergeValue, "UTF-8");


                // Default for merge is "false" for session files,  "true" otherwise
                boolean merge;
                if (mergeValue != null) {
                    // Explicit setting
                    merge = mergeValue.equalsIgnoreCase("true");
                } else if (file.endsWith(".xml") || file.endsWith(".php") || file.endsWith(".php3")) {
                    // Session file
                    merge = false;
                } else {
                    // Data file
                    merge = true;
                }

                String name = params.get("name");
                String format = params.get("format");
                String locus = params.get("locus");
                String index = params.get("index");
                String coverage = params.get("coverage");
                String sort = params.get("sort");
                String sortTag = params.get("sortTag");
                result = cmdExe.loadFiles(file, index, coverage, name, format, locus, merge, params, sort, sortTag);
            } else {
                return ("ERROR Parameter \"file\" is required");
            }
        } else if (command.equals("/reload") || command.equals("/goto")) {
            String locus = params.get("locus");
            IGV.getInstance().goToLocus(locus);
        } else if (command.equals("/execute")) {
            String param = StringUtils.decodeURL(params.get("command"));
            return cmdExe.execute(param);
        } else {
            return ("ERROR Unknown command: " + command);
        }

        return result;
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

                //This might look backwards, but it isn't.
                //Parameters must be URL encoded, including the file parameter
                //CommandExecutor URL-decodes the file parameter sometimes, but not always
                //So we URL-decode iff CommandExecutor doesn't
                boolean cmdExeWillDecode = (fileParams.contains(key) || indexParams.contains(key)) && CommandExecutor.needsDecode(kv[1]);

                String value = cmdExeWillDecode ? kv[1] : StringUtils.decodeURL(kv[1]);
                params.put(kv[0], value);
            }
        }
        return params;

    }


    /**
     * Compute a socket key according to the WebSocket RFC.  This method is here because this is the only class that uses it.
     *
     * @param input
     * @return
     * @throws NoSuchAlgorithmException
     */
    static String computeResponseKey(String input) throws NoSuchAlgorithmException, UnsupportedEncodingException {

        java.security.MessageDigest digest = null;

        digest = java.security.MessageDigest.getInstance("SHA-1");

        digest.reset();

        digest.update(input.getBytes("UTF-8"));

        return new String(Base64Coder.encode(digest.digest()));
    }
}
