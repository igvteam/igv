/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.action;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.logging.*;
import org.broad.igv.session.Session;
import org.broad.igv.session.SessionWriter;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.LongRunningTask;

import java.awt.event.ActionEvent;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * @author jrobinso
 */
public class ReloadTracksMenuAction extends MenuAction {

    static Logger log = LogManager.getLogger(SaveSessionMenuAction.class);
    IGV igv;

    /**
     * @param label
     * @param mnemonic
     * @param igv
     */
    public ReloadTracksMenuAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
    }

    /**
     * Method description
     *
     * @param e
     */
    @Override
    public void actionPerformed(ActionEvent e) {

        String currentSessionFilePath = igv.getSession().getPath();
        Session currentSession = igv.getSession();
        currentSession.setPath(currentSessionFilePath);
        String xml = (new SessionWriter()).createXmlFromSession(currentSession, null);

        igv.resetSession(currentSessionFilePath);
        final InputStream inputStream = new ByteArrayInputStream(xml.getBytes());

        Runnable runnable = () -> {
            try {
                igv.loadSessionFromStream(currentSessionFilePath, inputStream);
            } catch (IOException ex) {
                log.error("Error reloading tracks", ex);
            }
        };
        LongRunningTask.submit(runnable);
    }
}
