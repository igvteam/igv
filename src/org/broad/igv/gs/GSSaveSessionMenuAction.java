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

package org.broad.igv.gs;

/**
 * @author Jim Robinson
 * @date 8/10/11
 */
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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

import org.apache.log4j.Logger;
import org.broad.igv.gs.GSFileBrowser;
import org.broad.igv.gs.dm.DMUtils;
import org.broad.igv.session.Session;
import org.broad.igv.session.SessionWriter;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.ui.action.MenuAction;
import org.broad.igv.ui.action.SaveSessionMenuAction;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.FileOutputStream;

/**
 * @author jrobinso
 */
public class GSSaveSessionMenuAction extends MenuAction {

    static Logger log = Logger.getLogger(SaveSessionMenuAction.class);

    // TODO -- The batch referenceFrame is likely to be used by many actions. Move this
    // member to a base class ?
    IGV mainFrame;

    /**
     * Constructs ...
     *
     * @param label
     * @param mainFrame
     */
    public GSSaveSessionMenuAction(String label, IGV mainFrame) {
        super(label);
        this.mainFrame = mainFrame;
        setToolTipText("Save session to GenomeSpace");
    }

    /**
     * Method description
     *
     * @param e
     */
    @Override
    public void actionPerformed(ActionEvent e) {

        WaitCursorManager.CursorToken token = null;

        try {

            String currentSessionFilePath = mainFrame.getSession().getPath();
            String initFile = currentSessionFilePath == null ? UIConstants.DEFAULT_SESSION_FILE : currentSessionFilePath;
            String sessionName = (new File(initFile)).getName();

            GSFileBrowser gsFileBrowser = new GSFileBrowser(IGV.getMainFrame(), GSFileBrowser.Mode.SAVE);
            gsFileBrowser.setVisible(true);
            String gsPath = gsFileBrowser.getPath();
            if (gsPath == null) return;
            if (!gsPath.endsWith(".xml")) {
                gsPath += ".xml";
            }

            File tmpDir = new File(System.getProperty("java.io.tmpdir"), System.getProperty("user.name"));
            if(!tmpDir.exists()) {
                tmpDir.mkdirs();
            }
            sessionName = (new File(gsPath)).getName();
            final File sessionFile = new File(tmpDir, sessionName);

            mainFrame.setStatusBarMessage("Saving session to " + sessionFile.getAbsolutePath());


            FileOutputStream out = null;
            token = WaitCursorManager.showWaitCursor();

            Session currentSession = mainFrame.getSession();
            currentSession.setPath(sessionFile.getAbsolutePath());
            (new SessionWriter()).saveSession(currentSession, sessionFile);

            DMUtils.uploadFile(sessionFile, gsPath);

            sessionFile.delete();


        } catch (Exception e2) {
            JOptionPane.showMessageDialog(mainFrame.getMainFrame(), "Error saving session file " + "(" + e2.getMessage() + ")");
            log.error("Failed to save session!", e2);
        } finally {
            if (token != null) WaitCursorManager.removeWaitCursor(token);
            mainFrame.resetStatusMessage();

        }


    }
}
