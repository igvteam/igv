package org.igv.jbrowse;

import org.igv.logging.*;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.Socket;
import java.net.UnknownHostException;

class CircviewSocketWriter {

    static Logger log = LogManager.getLogger(CircviewSocketWriter.class);

    static String send(String json) {
        return send(json, false);
    }

    static String send(String json, boolean suppressErrors) {
        Socket socket = null;
        PrintWriter out = null;
        BufferedReader in = null;
        try {
            String host = PreferencesManager.getPreferences().get(Constants.CIRC_VIEW_HOST);
            int port = PreferencesManager.getPreferences().getAsInt(Constants.CIRC_VIEW_PORT);
            socket = new Socket(host, port);
            out = new PrintWriter(socket.getOutputStream(), true);
            in = new BufferedReader(new InputStreamReader(socket.getInputStream()));

            out.println(json);
            out.flush();
            String response = in.readLine();
            return response;
        } catch (UnknownHostException e) {
            String err = "Unknown host exception: " + e.getMessage();
            if (!suppressErrors) {
                log.error(e);
            }
            return err;

        } catch (IOException e) {
            String message = "IO Exception: " + e.getMessage();
            if (!suppressErrors) {
                log.error(message, e);
            }
            return message;
        } finally {
            try {
                in.close();
                out.close();
                socket.close();
            } catch (IOException e) {
                log.error(e);
            }
        }
    }
}
