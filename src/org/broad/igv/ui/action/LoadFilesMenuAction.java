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
package org.broad.igv.ui.action;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.util.FileChooserDialog;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 */
public class LoadFilesMenuAction extends MenuAction {

    static Logger log = Logger.getLogger(LoadFilesMenuAction.class);
    IGV mainFrame;

    public LoadFilesMenuAction(String label, int mnemonic, IGV mainFrame) {
        super(label, null, mnemonic);
        this.mainFrame = mainFrame;
        setToolTipText(UIConstants.LOAD_TRACKS_TOOLTIP);
    }

    @Override
    public void actionPerformed(ActionEvent e) {

        UIUtilities.invokeOnEventThread(new Runnable() {

            public void run() {
                loadFiles(chooseTrackFiles());
            }
        });
    }

    private File[] chooseTrackFiles() {

        File lastDirectoryFile =
                PreferenceManager.getInstance().getLastTrackDirectory();

        // Get Track Files
        FileChooserDialog trackFileDialog = mainFrame.getTrackFileChooser();
        trackFileDialog.setLocationRelativeTo(mainFrame.getMainFrame());
        trackFileDialog.setTitle("Select Files");
        trackFileDialog.setMultiSelectionEnabled(true);
        trackFileDialog.setSelectedFile(null);
        trackFileDialog.setCurrentDirectory(lastDirectoryFile);
        trackFileDialog.setFileSelectionMode(JFileChooser.FILES_ONLY);


        trackFileDialog.setVisible(true);

        // Tracks.  Simulates multi-file select
        File[] trackFiles = null;

        if (!trackFileDialog.isCanceled()) {

            File lastFile = trackFileDialog.getSelectedFile();
            if (lastFile != null && !lastFile.isDirectory()) {
                lastFile = lastFile.getParentFile();

                // Store the last accessed file location
                PreferenceManager.getInstance().setLastTrackDirectory(lastFile);
            }

            trackFiles = trackFileDialog.getSelectedFiles();
            trackFileDialog.setSelectedFile(null);
        }

        mainFrame.resetStatusMessage();

        return trackFiles;
    }

    private void loadFiles(File[] files) {

        if (files != null && files.length > 0) {

            List<File> validFileList = new ArrayList();
            StringBuffer buffer = new StringBuffer();
            buffer.append("File(s) not found: ");
            boolean allFilesExist = true;
            for (File file : files) {

                if (!file.exists()) {
                    allFilesExist = false;
                    buffer.append("\n\t");
                    buffer.append(file.getAbsolutePath());
                } else {

                    String path = file.getAbsolutePath();
                    if (path.endsWith(Globals.SESSION_FILE_EXTENSION)) {
                        // TODO -- a better test for session file than just the extension!
                        MessageUtils.showAndLogErrorMessage(mainFrame.getMainFrame(),
                                "File " + path +
                                        " appears to be an IGV Session file - " +
                                        "please use the Open Session menu item " +
                                        "to load it.", log);
                    } else {
                        validFileList.add(file);
                    }
                }

            }
            files = validFileList.toArray(new File[validFileList.size()]);

            if (!allFilesExist) {
                MessageUtils.showAndLogErrorMessage(mainFrame.getMainFrame(), buffer.toString(), log);
            }

            if (files.length > 0) {

                // Create DataResouceLocators for the selected files
                final List<ResourceLocator> locators = new ArrayList(files.length);
                for (File f : files) {
                    locators.add(new ResourceLocator(f.getAbsolutePath()));
                }

                mainFrame.loadTracks(locators);
            }
        }
    }
}
