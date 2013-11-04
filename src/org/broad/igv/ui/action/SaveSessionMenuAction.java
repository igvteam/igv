/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.action;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.session.Session;
import org.broad.igv.session.SessionWriter;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.ui.util.FileDialogUtils;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.IOException;

/**
 * @author jrobinso
 */
public class SaveSessionMenuAction extends MenuAction {

    static Logger log = Logger.getLogger(SaveSessionMenuAction.class);
    IGV igv;

    /**
     *
     *
     * @param label
     * @param mnemonic
     * @param igv
     */
    public SaveSessionMenuAction(String label, int mnemonic, IGV igv) {
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


        File sessionFile = null;

        String currentSessionFilePath = igv.getSession().getPath();

        String initFile = currentSessionFilePath == null ? UIConstants.DEFAULT_SESSION_FILE : currentSessionFilePath;
        sessionFile = FileDialogUtils.chooseFile("Save Session",
                PreferenceManager.getInstance().getLastSessionDirectory(),
                new File(initFile),
                FileDialogUtils.SAVE);


        if (sessionFile == null) {
            igv.resetStatusMessage();
            return;
        }


        String filePath = sessionFile.getAbsolutePath();
        if (!filePath.toLowerCase().endsWith(".xml")) {
            sessionFile = new File(filePath + ".xml");
        }

        igv.setStatusBarMessage("Saving session to " + sessionFile.getAbsolutePath());


        final File sf = sessionFile;
        WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();
        try {
            saveSession(igv, sf);
            // No errors so save last location
            PreferenceManager.getInstance().setLastSessionDirectory(sf.getParentFile());

        } catch (Exception e2) {
            JOptionPane.showMessageDialog(igv.getMainFrame(), "There was an error writing to " + sf.getName() + "(" + e2.getMessage() + ")");
            log.error("Failed to save session!", e2);
        } finally {
            WaitCursorManager.removeWaitCursor(token);
            igv.resetStatusMessage();


        }


    }

    /**
     * Saves current IGV session to {@code targetFile}. As a side effect,
     * sets the current sessions path (does NOT set the last session directory)
     * @param igv
     * @param targetFile
     * @throws IOException
     */
    public static void saveSession(IGV igv, File targetFile) throws IOException{
        Session currentSession = igv.getSession();
        currentSession.setPath(targetFile.getAbsolutePath());
        (new SessionWriter()).saveSession(currentSession, targetFile);
    }
}
