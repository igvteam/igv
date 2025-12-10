/*
This code was specifically designed to work on Unix-like systems
where the 'ls' command is available. It extends the JFileChooser FileDialog
It was written to handle directories with huge numbers of subdirectories/files
more efficiently by leveraging the native 'ls' command for listing contents.

Note this only used if environment variable USE_CUSTOM_FILEDIALOG is set to "true".

Background: We are running IGV in an XPRA desktop session on a remote server.
The filesystem has a directory mounted where the backend is a cloud object store (s3)
Some of the folders have more than 3000 files/folders in them and the standard FielDialog either took
hours to fetch the list of files/folders. At the same time we observed that running 'ls' on these folders
was nearly instantaneous (a few seconds at most).
*/
package org.broad.igv.ui.util;

import org.broad.igv.logging.*;

import javax.swing.*;
import java.io.File;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.InputStreamReader;

public class LsFileChooser extends JFileChooser {

    private static Logger log = LogManager.getLogger(LsFileChooser.class);

    public LsFileChooser(File currentDirectory) {
        super(currentDirectory);

        // Listen for directory changes
        this.addPropertyChangeListener(JFileChooser.DIRECTORY_CHANGED_PROPERTY, evt -> {
            File dir = getCurrentDirectory();
            if (dir != null && dir.isDirectory()) {
                scanDirectory(dir);
            }
        });
    }

    @Override
    // This method can be called to refresh the current directory using the same logic as when the directory is changed
    public void rescanCurrentDirectory() {
        File dir = getCurrentDirectory();
        if (dir != null && dir.isDirectory()) {
            scanDirectory(dir);
        }
    }

    /**
     * Scans the given directory and logs the number of files and directories.
     *
     * @param dir The directory to scan.
     */
    private void scanDirectory(File dir) {
        int fileCount = 0;
        int dirCount = 0;
        try {
            Process process = new ProcessBuilder("ls", "-1", dir.getAbsolutePath()).start();
            try (BufferedReader reader = new BufferedReader(
                    new InputStreamReader(process.getInputStream()))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    File file = new File(dir, line);
                    if (file.isDirectory()) {
                        dirCount++;
                    } else {
                        fileCount++;
                    }
                }
            }
            int exitCode = process.waitFor();
            if (exitCode != 0) {
                log.error("Error: ls command failed with exit code " + exitCode);
                return;
            }
            log.info("Directory: " + dir.getAbsolutePath() +
                    " | Files: " + fileCount + " | Folders: " + dirCount);
        } catch (IOException | InterruptedException e) {
            log.error("Error scanning directory: " + dir.getAbsolutePath(), e);
        }
    }
}
