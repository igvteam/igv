/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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
package org.broad.igv.main;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.IGVMainFrame;
import org.broad.igv.ui.WaitCursorManager;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.ServerSocket;
import java.net.Socket;
import java.net.URLDecoder;
import java.util.HashMap;
import java.util.Map;

public class CommandListener implements Runnable {

    private static Logger log = Logger.getLogger(CommandListener.class);

    private static CommandListener listener;

    private int port = -1;

    private CommandListener(int port) {
        this.port = port;
    }

    public static synchronized void halt() {
        if (listener != null) {
            listener.halt = true;
        }
    }

    public static synchronized void start(int port) {
        listener = new CommandListener(port);
        Thread listenerThread = new Thread(listener);
        listenerThread.start();
    }

    boolean halt = false;

    public void run() {

        CommandExecutor cmdExe = new CommandExecutor();

        ServerSocket serverSocket = null;
        try {
            serverSocket = new ServerSocket(port);
            log.info("Listening on port " + port);

            while (halt == false) {
                Socket clientSocket = null;
                try {
//                    System.out.println("Listening on port: " + port);
                    clientSocket = serverSocket.accept();

                    PrintWriter out = new PrintWriter(clientSocket.getOutputStream(), true);
                    BufferedReader in = new BufferedReader(new InputStreamReader(clientSocket.getInputStream()));

                    String inputLine;


                    while ((inputLine = in.readLine()) != null) {
                        Globals.batch = true;
                        Globals.setSuppressMessages(true);
                        WaitCursorManager.CursorToken cursorToken = null;
                        try {
                            String cmd = inputLine;
                            if (halt == true) {
                                if (cmd.startsWith("GET")) {
                                    sendHTTPResponse(out, "ERROR IGV port is closed");
                                } else {
                                    out.println("ERROR IGV port is closed");
                                }
                                break;
                            }
                            if (cmd.startsWith("GET")) {
                                String result = processGet(cmd, in, cmdExe);
                                sendHTTPResponse(out, result);
                                clientSocket.close();

                                clientSocket = serverSocket.accept();
                                out = new PrintWriter(clientSocket.getOutputStream(), true);
                                in = new BufferedReader(new InputStreamReader(clientSocket.getInputStream()));

                            } else {
                                out.println(cmdExe.execute(inputLine));
                            }
                        } finally {
                            Globals.setSuppressMessages(false);
                            Globals.batch = false;
                            if (cursorToken != null) WaitCursorManager.removeWaitCursor(cursorToken);
                        }

                    }


                    out.close();
                    in.close();
                    clientSocket.close();
                } catch (IOException e) {
                    log.error("Accept failed.", e);
                } finally {
                    clientSocket.close();
                    Globals.setSuppressMessages(false);
                }

            }

        } catch (java.net.BindException e) {
            log.error(e);
        } catch (Exception e) {
            log.error("Could not listen on port: " + port, e);
        } finally {
            try {
                if (serverSocket != null) {
                    serverSocket.close();
                }
            } catch (IOException ex) {
                log.error("Error closing command listener socket", ex);
            }
        }
    }

    private void sendHTTPResponse(PrintWriter out, String result) {
        out.println("HTTP/1.0 204 OK");
        out.println(" Server: IGV");
        out.println("Connection: close");
        out.println();
        out.println(result);
        out.println();
        out.close();
    }

    /**
     * Process an http get request.
     *
     * @param line
     * @param reader
     * @return
     * @throws IOException
     */

    private String processGet(String line, BufferedReader reader, CommandExecutor cmdExe) throws IOException {

        String nextLine = URLDecoder.decode(line);
        String result = "OK";

        String[] tokens = nextLine.split(" ");
        if (tokens.length < 2) {
            return "ERROR unexpected command line: " + line;
        } else {
            String[] parts = tokens[1].split("\\?");
            if (parts.length < 2) {
                return ("ERROR unexpected command line: " + line);
            } else {
                String command = parts[0];
                Map<String, String> params = parseParameters(parts[1]);
                final IGVMainFrame mainFrame = IGVMainFrame.getInstance();

                // Trick to force window to front, the setAlwaysOnTop works on a Mac,  toFront() does nothing.
                mainFrame.toFront();
                mainFrame.setAlwaysOnTop(true);
                mainFrame.setAlwaysOnTop(false);

                if (command.equals("/load")) {
                    if (params.containsKey("file")) {
                        String genomeID = params.get("genome");
                        String mergeValue = params.get("merge");
                        String locus = params.get("locus");
                        if (genomeID != null) {
                            mainFrame.selectGenomeFromList(genomeID);
                        }

                        // Default for merge is "true"
                        boolean merge = mergeValue == null ? true : mergeValue.equalsIgnoreCase("true");

                        result = cmdExe.execute("hget " + params.get("file") + " " + locus + " " + merge);
                    } else {
                        return ("ERROR Parameter \"file\" is required");
                    }
                } else if (command.equals("/reload") || command.equals("/goto")) {
                    String locus = params.get("locus");
                    mainFrame.goToLocus(locus);
                } else {
                    return ("ERROR Unknown command: " + command);
                }
            }
        }

        // Consume the remainder of the request, if any.  This is important to free the connection.
        while (nextLine.length() > 0) {
            nextLine = reader.readLine();
        }
        return result;
    }

    private Map<String, String> parseParameters(String parameterString) {
        HashMap<String, String> params = new HashMap();
        String[] kvPairs = parameterString.split("&");
        for (String kvString : kvPairs) {
            String[] kv = kvString.split("=");
            if (kv.length == 1) {
                params.put(kv[0], null);
            } else {
                params.put(kv[0], kv[1]);
            }
        }
        return params;

    }
}
