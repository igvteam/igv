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

package org.broad.igv.gs;

/**
 * @author Jim Robinson
 * @date 8/10/11
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
