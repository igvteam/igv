package org.igv.ui.action;


import org.igv.logging.*;
import org.igv.session.Session;
import org.igv.session.SessionWriter;
import org.igv.ui.IGV;
import org.igv.util.LongRunningTask;

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
