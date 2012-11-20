/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
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
package org.broad.igv.gs;

import org.apache.log4j.Logger;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.action.MenuAction;
import org.broad.igv.ui.util.MessageUtils;

import javax.swing.*;
import java.awt.event.ActionEvent;

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
            if (IGV.getInstance().getAllTracks().size() > 0) {
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
                        boolean merge = false;
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

