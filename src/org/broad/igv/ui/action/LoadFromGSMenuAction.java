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

package org.broad.igv.ui.action;

import org.apache.log4j.Logger;
import org.broad.igv.gs.GSFileBrowser;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ResourceLocator;

import java.awt.event.ActionEvent;
import java.util.Arrays;

/**
 * @author jrobinso
 * @date Jun 5, 2011
 */
public class LoadFromGSMenuAction extends MenuAction {

    static Logger log = Logger.getLogger(LoadFilesMenuAction.class);
    private IGV mainFrame;

    public LoadFromGSMenuAction(String label, int mnemonic, IGV mainFrame) {
        super(label, null, mnemonic);
        this.mainFrame = mainFrame;
        setToolTipText("Load from GenomeSpace");
    }

    @Override
    public void actionPerformed(ActionEvent e) {


        try {
            GSFileBrowser dlg = new GSFileBrowser(IGV.getMainFrame());
            dlg.setVisible(true);

            String url = dlg.getFileURL();
            if (url != null) {
                IGV.getInstance().loadTracks(Arrays.asList(new ResourceLocator(url)));
            }
        } catch (Exception e1) {
            log.error("Error fetching directory listing on GenomeSpace server.", e1);
            MessageUtils.showMessage("Error fetching directory listing on GenomeSpace server: " + e1.getMessage());
        }

        /*JPanel ta = new JPanel();
        ta.setPreferredSize(new Dimension(600, 20));
        if (e.getActionCommand().equalsIgnoreCase(LOAD_FROM_URL)) {
            String url = JOptionPane.showInputDialog(IGV.getMainFrame(), ta, "Enter URL (http or ftp)", JOptionPane.QUESTION_MESSAGE);
            if (url != null && url.trim().length() > 0) {
                if (url.endsWith(".xml")) {
                    try {
                        mainFrame.doRestoreSession(new URL(url), null);
                    } catch (Exception ex) {
                        MessageUtils.showMessage("Error loading url: " + url + " (" + ex.toString() + ")");
                    }
                } else {
                    ResourceLocator rl = new ResourceLocator(url.trim());
                    mainFrame.loadTracks(Arrays.asList(rl));

                }
            }
        } else if ((e.getActionCommand().equalsIgnoreCase(LOAD_FROM_DAS))) {
            String url = JOptionPane.showInputDialog(IGV.getMainFrame(), ta, "Enter DAS feature source URL",
                    JOptionPane.QUESTION_MESSAGE);
            if (url != null && url.trim().length() > 0) {
                ResourceLocator rl = new ResourceLocator(url.trim());
                rl.setType("das");
                mainFrame.loadTracks(Arrays.asList(rl));
            }
        }
        */
    }

}

