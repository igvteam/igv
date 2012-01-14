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
package org.broad.igv.gs;

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.action.MenuAction;
import org.broad.igv.ui.util.FileChooserDialog;
import org.broad.igv.ui.util.MessageUtils;
import org.json.JSONException;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.IOException;
import java.net.URL;

/**
 * This menu action classes is used for both "Open Session ..." and load recent
 * session menu items.  In the "Open Session..." the user has to specify a
 * session file through the file menu.  For "load recent"  the action is
 * instantiated with a specific session file.
 *
 * @author jrobinso
 */
public class GSOpenSessionMenuAction extends MenuAction {

    private static Logger log = Logger.getLogger(GSOpenSessionMenuAction.class);
    private IGV mainFrame;
    private String sessionFile = null;
    private boolean autoload = false;

    public GSOpenSessionMenuAction(String label, IGV mainFrame) {
        super(label);
        this.mainFrame = mainFrame;
        setToolTipText(UIConstants.RESTORE_SESSION_TOOLTIP);
    }

    @Override
    public void actionPerformed(ActionEvent e) {

        try {
            GSFileBrowser fileBrowser = new GSFileBrowser(IGV.getMainFrame(), GSFileBrowser.Mode.OPEN);
            fileBrowser.setVisible(true);

            sessionFile = fileBrowser.getFileURL();
            if (sessionFile == null || sessionFile.length() == 0) {
                return;
            }

            // If anything has been loaded warn the users.  Popping up the
            // warning all the time will get annoying.
            if (IGV.getInstance().getAllTracks(false).size() > 0) {
                int status =
                        JOptionPane.showConfirmDialog(
                                mainFrame.getMainFrame(),
                                UIConstants.OVERWRITE_SESSION_MESSAGE,
                                null,
                                JOptionPane.OK_CANCEL_OPTION,
                                JOptionPane.PLAIN_MESSAGE,
                                null);

                if (status == JOptionPane.CANCEL_OPTION ||
                        status == JOptionPane.CLOSED_OPTION) {
                    return;
                }
            }

            if (sessionFile != null) {

                SwingWorker worker = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        boolean merge=false;
                        mainFrame.doRestoreSession(sessionFile, null, merge);
                        return null;
                    }
                };
                worker.execute();

            }
        } catch (Exception e1) {
            log.error("Error restoring session", e1);
            MessageUtils.confirm("Error loading session from GenomeSpace.  " + e1.toString());
        }
    }

}

