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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.action;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.session.Session;
import org.broad.igv.session.SessionWriter;
import org.broad.igv.ui.IGVMainFrame;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.ui.util.FileDialogUtils;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * @author jrobinso
 */
public class SaveSessionMenuAction extends MenuAction {

    static Logger log = Logger.getLogger(SaveSessionMenuAction.class);

    // TODO -- The main referenceFrame is likely to be used by many actions. Move this
    // member to a base class ?
    IGVMainFrame mainFrame;

    /**
     * Constructs ...
     *
     * @param label
     * @param mnemonic
     * @param mainFrame
     */
    public SaveSessionMenuAction(String label, int mnemonic, IGVMainFrame mainFrame) {
        super(label, null, mnemonic);
        this.mainFrame = mainFrame;
        setToolTipText(UIConstants.SAVE_SESSION_TOOLTIP);
    }

    /**
     * Method description
     *
     * @param e
     */
    @Override
    public void actionPerformed(ActionEvent e) {


        File sessionFile = null;

        String currentSessionFilePath = mainFrame.getSession().getPath();

        // If no previous session file or we're doing a save as
        if (true) // currentSessionFilePath == null)
        {

            String initFile = currentSessionFilePath == null ? UIConstants.DEFAULT_SESSION_FILE : currentSessionFilePath;
            sessionFile = FileDialogUtils.chooseFile("Save Session",
                    PreferenceManager.getInstance().getLastSessionDirectory(),
                    new File(initFile),
                    FileDialogUtils.SAVE);


            if (sessionFile == null) {
                mainFrame.resetStatusMessage();
                return;
            }


            String filePath = sessionFile.getAbsolutePath();
            if (!filePath.toLowerCase().endsWith(Globals.SESSION_FILE_EXTENSION)) {
                sessionFile = new File(filePath + Globals.SESSION_FILE_EXTENSION);
            }

        } else {
            sessionFile = new File(currentSessionFilePath);
        }

        mainFrame.setStatusBarMessage("Saving session to " + sessionFile.getAbsolutePath());

        // Get rid of old file before creating the new one
        if (sessionFile.exists()) {
            sessionFile.delete();
        }

        final File sf = sessionFile;

        FileOutputStream out = null;
        WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();
        try {

            Session currentSession = mainFrame.getSession();
            currentSession.setFilePath(sf.getAbsolutePath());
            (new SessionWriter()).saveSession(currentSession, sf);

            // No errors so save last location
            PreferenceManager.getInstance().setLastSessionDirectory(sf.getParentFile());

        } catch (Exception e2) {
            JOptionPane.showMessageDialog(mainFrame, "There was an error writing to " + sf.getName() + "(" + e2.getMessage() + ")");
            log.error("Failed to save session!", e2);
        } finally {
            WaitCursorManager.removeWaitCursor(token);
            mainFrame.resetStatusMessage();

            if (out != null) {

                try {
                    out.close();
                } catch (IOException exception) {
                    log.error("Failed to close session file!", exception);
                }
            }
        }


    }
}
