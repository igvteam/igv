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
