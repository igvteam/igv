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
import org.broad.igv.ui.IGV;

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
     * Loop forever, processing client requests synchronously.  The server is single threaded, because in most cases
     * we would not know how to process commands ssychronously
     */
    public void run() {

        CommandExecutor cmdExe = new CommandExecutor();

        try {
            serverSocket = new ServerSocket(port);
            log.info("Listening on port " + port);

            while (true) {
                clientSocket = serverSocket.accept();
                processClientSession(cmdExe);
                if (clientSocket != null) {
                    try {
                        clientSocket.close();
                        clientSocket = null;
                    } catch (IOException e) {
                        e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                    }
                }
            }


        } catch (java.net.BindException e) {
            log.error(e);
        } catch (ClosedByInterruptException e) {
            // Not really an error, caused by an interrupt

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
                    String result = processGet(cmd, in, cmdExe);
                    sendHTTPResponse(out, result);
                } else {
                    Globals.setBatch(true);
                    Globals.setSuppressMessages(true);
                    out.println(cmdExe.execute(inputLine));
                }
            }
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
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
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }

        if (serverSocket != null) {
            try {
                serverSocket.close();
                serverSocket = null;
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
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

        String nextLine = line;
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
                final Frame mainFrame = IGV.getMainFrame();

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
                            IGV.getFirstInstance().selectGenomeFromList(genomeID);
                        }

                        // Default for merge is "true"
                        boolean merge = mergeValue == null ? true : mergeValue.equalsIgnoreCase("true");

                        result = cmdExe.execute("hget " + params.get("file") + " " + locus + " " + merge);
                    } else {
                        return ("ERROR Parameter \"file\" is required");
                    }
                } else if (command.equals("/reload") || command.equals("/goto")) {
                    String locus = params.get("locus");
                    IGV.getFirstInstance().goToLocus(locus);
                } else {
                    return ("ERROR Unknown command: " + command);
                }
            }
        }

        // Consume the remainder of the request, if any.  This is important to free the connection.
        while (nextLine != null && nextLine.length() > 0) {
            nextLine = reader.readLine();
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
        HashMap<String, String> params = new HashMap();
        String[] kvPairs = parameterString.split("&");
        for (String kvString : kvPairs) {
            String[] kv = kvString.split("=");
            if (kv.length == 1) {
                params.put(kv[0], null);
            } else {
                String key = kv[0];
                // Special treatment of locus string, need to preserve encoding of spaces
                String value = key.equals("locus") ? kv[1] : URLDecoder.decode(kv[1]);
                params.put(kv[0], value);
            }
        }
        return params;

    }
}
