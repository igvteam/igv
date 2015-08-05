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
package org.broad.igv.ui.action;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ResourceLocator;

import java.awt.event.ActionEvent;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 */
public class LoadFilesMenuAction extends MenuAction {

    static Logger log = Logger.getLogger(LoadFilesMenuAction.class);
    IGV igv;

    public LoadFilesMenuAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
    }

    @Override
    public void actionPerformed(ActionEvent e) {

                loadFiles(chooseTrackFiles());

    }

    private File[] chooseTrackFiles() {

        File lastDirectoryFile = PreferenceManager.getInstance().getLastTrackDirectory();

        // Get Track Files
        final PreferenceManager prefs = PreferenceManager.getInstance();

        // Tracks.  Simulates multi-file select
        File[] trackFiles = FileDialogUtils.chooseMultiple("Select Files", lastDirectoryFile, null);

        if (trackFiles != null && trackFiles.length > 0) {

            File lastFile = trackFiles[0];
            if (lastFile != null) {
                PreferenceManager.getInstance().setLastTrackDirectory(lastFile);
            }
        }
        igv.resetStatusMessage();
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
                        final String msg = "File " + path +
                                " appears to be an IGV Session file - " +
                                "please use the Open Session menu item " +
                                "to load it.";
                        log.error(msg);
                        MessageUtils.showMessage(msg);
                    } else {
                        validFileList.add(file);
                    }
                }

            }
            files = validFileList.toArray(new File[validFileList.size()]);

            if (!allFilesExist) {
                final String msg = buffer.toString();
                log.error(msg);
                MessageUtils.showMessage(msg);
            }

            if (files.length > 0) {

                // Create DataResouceLocators for the selected files
                final List<ResourceLocator> locators = new ArrayList(files.length);
                for (File f : files) {
                    locators.add(new ResourceLocator(f.getAbsolutePath()));
                }

                igv.loadTracks(locators);
            }
        }
    }
}
