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

import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.session.SessionReader;
import org.broad.igv.track.AttributeManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.ResourceLocator;

import java.awt.event.ActionEvent;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author jrobinso
 */
public class LoadFilesMenuAction extends MenuAction {

    private static final Logger log = LogManager.getLogger(LoadFilesMenuAction.class);


    public enum Type {TRACK, SESSION, SAMPLE_INFO}

    private final IGV igv;
    private final Type type;

    public LoadFilesMenuAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
        this.type = Type.TRACK;
    }

    public LoadFilesMenuAction(String label, int mnemonic, IGV igv, Type type) {
        super(label, null, mnemonic);
        this.igv = igv;
        this.type = type;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        loadFiles(chooseTrackFiles());
    }

    private File[] chooseTrackFiles() {

        File lastDirectoryFile = PreferencesManager.getPreferences().getLastTrackDirectory();

        // Tracks.  Simulates multi-file select
        File[] trackFiles = FileDialogUtils.chooseMultiple("Select Files", lastDirectoryFile, null);

        if (trackFiles != null && trackFiles.length > 0) {

            File lastFile = trackFiles[0];
            if (lastFile != null) {
                PreferencesManager.getPreferences().setLastTrackDirectory(lastFile);
            }
        }
        igv.resetStatusMessage();
        return trackFiles;
    }

    private void loadFiles(final File[] files) {

        if (files != null && files.length > 0) {

            final List<File> validFiles = new ArrayList<>();
            final List<File> missingFiles = new ArrayList<>();

            for (File file : files) {
                if (!file.exists()) {
                    missingFiles.add(file);
                } else {
                    String path = file.getAbsolutePath();
                    if (this.type == Type.SAMPLE_INFO) {
                        LongRunningTask.submit(() -> {
                            AttributeManager.getInstance().loadSampleInfo(new ResourceLocator(path));
                            igv.revalidateTrackPanels();
                        });
                    } else if (SessionReader.isSessionFile(path)) {
                        LongRunningTask.submit(() -> this.igv.loadSession(path, null));
                    } else {
                        validFiles.add(file);
                    }
                }
            }

            if (!missingFiles.isEmpty()) {
                String msg = missingFiles.stream()
                        .map(File::getAbsolutePath)
                        .collect(Collectors.joining("\n\t", "File(s) not found: \n\t", ""));
                log.error(msg);
                MessageUtils.showMessage(msg);
            }

            if (!validFiles.isEmpty()) {
                // Create DataResourceLocators for the selected files
                final List<ResourceLocator> locators = ResourceLocator.getLocators(validFiles);
                igv.addToRecentUrls(locators);
                igv.loadTracks(locators);
            }
        }
    }
}
