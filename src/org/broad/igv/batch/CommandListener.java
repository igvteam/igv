/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */
package org.broad.igv.batch;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.StringUtils;

import java.awt.*;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.ServerSocket;
import java.net.Socket;
import java.net.URLDecoder;
import java.nio.channels.ClosedByInterruptException;
import java.util.HashMap;
import java.util.Map;
import org.broad.igv.PreferenceManager;

public class CommandListener implements Runnable {

    private static Logger log = Logger.getLogger(CommandListener.class);

    private static CommandListener listener;

    private int port = -1;
    private ServerSocket serverSocket = null;
    private Socket clientSocket = null;
    private Thread listenerThread;
    boolean halt = false;

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
                    String command = null;
                    Map<String, String> params = null;
                    String[] tokens = inputLine.split(" ");
                    if (tokens.length < 2) {
                        sendHTTPResponse(out, "ERROR unexpected command line: " + inputLine);
                        return;
                    } else {
                        String[] parts = tokens[1].split("\\?");
                        if (parts.length < 2) {
                            sendHTTPResponse(out, "ERROR unexpected command line: " + inputLine);
                            return;
                        } else {
                            command = parts[0];
                            params = parseParameters(parts[1]);
                        }
                    }

                    // Consume the remainder of the request, if any.  This is important to free the connection.
                    String nextLine = in.readLine();
                    while (nextLine != null && nextLine.length() > 0) {
                        nextLine = in.readLine();
                    }

                    // If a callback (javascript) function is specified write it back immediately.  This function
                    // is used to cancel a timeout handler
                    String callback = params.get("callback");
                    if (callback != null) {
                        sendHTTPResponse(out, callback);
                    }

                    processGet(command, params, cmdExe);

                    // If no callback was specified write back a "no response" header
                    if (callback == null) {
                        sendHTTPResponse(out, null);
                    }

                    // http sockets are used for one request onle
                    return;

                } else {
                    // Port command
                    Globals.setBatch(true);
                    Globals.setSuppressMessages(true);
                    final String response = cmdExe.execute(inputLine);
                    out.println(response);
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

    private static final String CRNL = "\r\n";
    private static final String CONTENT_TYPE = "Content-Type: ";
    private static final String HTTP_RESPONSE = "HTTP/1.1 200 OK";
    private static final String HTTP_NO_RESPONSE = "HTTP/1.1 204 No Response";
    private static final String CONTENT_LENGTH = "Content-Length: ";
    private static final String CONTENT_TYPE_TEXT_HTML = "text/html";
    private static final String CONNECTION_CLOSE = "Connection: close";
    private static final String NO_CACHE = "Cache-Control: no-cache, no-store";

    private void sendHTTPResponse(PrintWriter out, String result) {

        out.println(result == null ? HTTP_NO_RESPONSE : HTTP_RESPONSE);
        if (result != null) {
            out.print(CONTENT_TYPE + CONTENT_TYPE_TEXT_HTML);
            out.print(CRNL);
            out.print(CONTENT_LENGTH + (result.length()));
            out.print(CRNL);
            out.print(NO_CACHE);
            out.print(CRNL);
            out.print(CONNECTION_CLOSE);
            out.print(CRNL);
            out.print(CRNL);
            out.print(result);
            out.print(CRNL);
        }
        out.close();
    }

    /**
     * Process an http get request.
     */

    private String processGet(String command, Map<String, String> params, CommandExecutor cmdExe) throws IOException {

        String result = "OK";
        final Frame mainFrame = IGV.getMainFrame();

        // Trick to force window to front, the setAlwaysOnTop works on a Mac,  toFront() does nothing.
        mainFrame.toFront();
        mainFrame.setAlwaysOnTop(true);
        mainFrame.setAlwaysOnTop(false);

        //track%20type=bigBed%20name=%27hES_HUES1_p28.RRBS_CpG_meth%27%20description=%27RRBS%20CpG%20methylation%20for%20hES_HUES1_p28.RRBS%27%20visibility=4%20useScore=1%20color=0,60,120

        /** from what server was IGV started? Used to link to other tools/apps */
        
        String server = params.get("server");
        if (server == null || server.trim().length()<1) server = PreferenceManager.getInstance().get(PreferenceManager.IONTORRENT_SERVER);
        else {
            PreferenceManager.getInstance().put(PreferenceManager.IONTORRENT_SERVER, server);
        }
        if (command.equals("/load")) {
            String file = params.get("file");
            if (file == null) {
                file = params.get("bigDataURL"); // <- UCSC track line
            }
            if (file == null) {
                file = params.get("sessionURL");  // <- older IGV option
            }
            if (file == null) {
                file = params.get("dataURL"); // <- Another UCSC option
            }


            if (file != null) {
                PreferenceManager.getInstance().put(PreferenceManager.IONTORRENT_RESULTS, file);
                String genome = params.get("genome");
                if (genome == null) {
                    genome = params.get("db");  // <- UCSC track line param
                }
                if (genome != null) {
                    genome = URLDecoder.decode(genome, "UTF-8");
                    if (IGV.getInstance().getSelectableGenomeIDs().contains(genome)) {
                        IGV.getInstance().selectGenomeFromList(genome);
                    }
                    else {
                       IGV.getInstance().loadGenome(genome.trim(), null);
                    }
                }

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

                String locus = params.get("locus");
                if (locus != null) {
                    locus = URLDecoder.decode(locus, "UTF-8");
                }
                result = cmdExe.loadFiles(file, locus, merge, name, params);
            } else {
                return ("ERROR Parameter \"file\" is required");
            }
        } else if (command.equals("/reload") || command.equals("/goto")) {
            String locus = params.get("locus");
            IGV.getInstance().goToLocus(locus);
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
                // Special treatment of locus string, need to preserve encoding of spaces
                String value = key.equals("locus") ? kv[1] : StringUtils.decodeURL(kv[1]);
                params.put(kv[0], value);
            }
        }
        return params;

    }
}
