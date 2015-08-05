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

