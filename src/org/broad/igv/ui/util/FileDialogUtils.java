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

package org.broad.igv.ui.util;

import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.filefilters.AlignmentFileFilter;
import org.broad.igv.ui.filefilters.CoverageFileFilter;

import javax.swing.*;
import javax.swing.filechooser.FileFilter;
import java.awt.*;
import java.io.File;
import java.io.FilenameFilter;

/**
 * @author jrobinso
 * @date Nov 17, 2010
 */
public class FileDialogUtils {

    public static int LOAD = FileDialog.LOAD;
    public static int SAVE = FileDialog.SAVE;


    public static File chooseFile(String title, File initialDirectory, int mode) {
        return chooseFile(title, initialDirectory, null, mode);
    }


    public static File chooseFile(String title) {
        return chooseFile(title, DirectoryManager.getUserDirectory(), null, FileDialog.LOAD);
    }

    public static File chooseFile(String title, File initialDirectory, File initialFile, int mode) {
        return chooseFile(title, initialDirectory, initialFile, null, JFileChooser.FILES_ONLY, mode);
    }

    public static File chooseFileOrDirectory(String title, File initialDirectory, File initialFile, int mode) {
        return chooseFile(title, initialDirectory, initialFile, null, JFileChooser.FILES_AND_DIRECTORIES, mode);
    }

    private static File chooseFile(String title, File initialDirectory, File initialFile, FilenameFilter filter,
                                   int directoriesMode, int mode) {

        if (initialDirectory == null && initialFile != null) {
            initialDirectory = initialFile.getParentFile();
        }

        // Strip off parent directory
        if (initialFile != null) initialFile = new File(initialFile.getName());

        // TODO -- use native dialogs for windows as well?
        if (Globals.IS_MAC && directoriesMode != JFileChooser.FILES_AND_DIRECTORIES) {
            return chooseNative(title, initialDirectory, initialFile, filter, directoriesMode, mode);
        } else {
            return chooseSwing(title, initialDirectory, initialFile, filter, directoriesMode, mode);
        }
    }

    public static File[] chooseMultiple(String title, File initialDirectory, final FilenameFilter filter) {
        JFileChooser fileChooser = getJFileChooser(title, initialDirectory, null, filter, JFileChooser.FILES_ONLY);
        fileChooser.setMultiSelectionEnabled(true);
        fileChooser.addChoosableFileFilter(new AlignmentFileFilter());
        fileChooser.addChoosableFileFilter(new CoverageFileFilter());
        // set the default file filter to "All"
        fileChooser.setFileFilter(fileChooser.getChoosableFileFilters()[0]);

        boolean approve = fileChooser.showOpenDialog(getParentFrame()) == JFileChooser.APPROVE_OPTION;
        if (approve) {
            return fileChooser.getSelectedFiles();
        } else {
            return null;
        }

    }

    public static File chooseDirectory(String title, File initialDirectory) {
        if (Globals.IS_MAC) {
            return chooseNative(title, initialDirectory, null, null, JFileChooser.DIRECTORIES_ONLY, LOAD);
        } else {
            return chooseSwing(title, initialDirectory, null, null, JFileChooser.DIRECTORIES_ONLY, LOAD);
        }
    }


    private static File chooseNative(String title, File initialDirectory, File initialFile, FilenameFilter filter,
                                     int directoryMode, int mode) {

        boolean directories = JFileChooser.DIRECTORIES_ONLY == directoryMode;
        System.setProperty("apple.awt.fileDialogForDirectories", String.valueOf(directories));
        Frame parentFrame = getParentFrame();
        FileDialog fd = new FileDialog(parentFrame, title);
        if (initialDirectory != null) {
            fd.setDirectory(initialDirectory.getAbsolutePath());
        }
        if (initialFile != null) {
            fd.setFile(initialFile.getName());
        }
        if (filter != null) {
            fd.setFilenameFilter(filter);
        }
        fd.setModal(true);
        fd.setMode(mode);
        fd.setVisible(true);

        String file = fd.getFile();
        String directory = fd.getDirectory();
        if (file != null && directory != null) {
            // Ugly MAC hack -- bug in their native file dialog
            if (Globals.IS_MAC && initialFile != null) {
                file = fixMacExtension(initialFile, file);
            }
            return new File(directory, file);
        } else {
            return null;
        }
    }

    private static File chooseSwing(String title, File initialDirectory, File initialFile, final FilenameFilter filter,
                                    int directoryMode, int mode) {

        UIManager.put("FileChooser.readOnly", Boolean.FALSE);
        JFileChooser fileChooser = getJFileChooser(title, initialDirectory, initialFile, filter, directoryMode);
        Frame parentFrame = getParentFrame();
        boolean approve;
        if (mode == LOAD) {
            approve = fileChooser.showOpenDialog(parentFrame) == JFileChooser.APPROVE_OPTION;
        } else {

            approve = fileChooser.showSaveDialog(parentFrame) == JFileChooser.APPROVE_OPTION;
        }

        if (approve) {
            return fileChooser.getSelectedFile();
        } else {
            return null;
        }

    }

    /**
     * @param title
     * @param initialDirectory
     * @param initialFile
     * @param filter
     * @param directoryMode    either JFileChooser.DIRECTORIES_ONLY, JFileChooser.FILES_ONLY, or
     *                         JFileChooser.DIRECTORIES_ONLY : JFileChooser.FILES_AND_DIRECTORIES
     * @return
     */
    private static JFileChooser getJFileChooser(String title, File initialDirectory, File initialFile,
                                                final FilenameFilter filter, int directoryMode) {
        JFileChooser fileChooser = new JFileChooser();
        if (initialDirectory != null) {
            fileChooser.setCurrentDirectory(initialDirectory);
        }
        if (initialFile != null) {
            fileChooser.setSelectedFile(initialFile);
        }
        if (filter != null) {
            fileChooser.setFileFilter(new FileFilter() {
                @Override
                public boolean accept(File file) {
                    return filter.accept(file.getParentFile(), file.getName());
                }

                @Override
                public String getDescription() {
                    return "";
                }
            });
        }

        fileChooser.setDialogTitle(title);
        fileChooser.setFileSelectionMode(directoryMode);

        return fileChooser;
    }


    /**
     * Fix for bug in MacOS "native" dialog.  If hide extension is on the extension is stripped from the dialog,
     * so far so good, but it is not added in the file returned.  So, if we know the extension expected from
     * initialFile add it back.  If not too bad.
     *
     * @param initialFile
     * @param fname
     */
    private static String fixMacExtension(File initialFile, String fname) {
        if (fname.contains(".")) {
            return fname;   // Has some sort of extension.  Should we compare to expected extension?
        }
        String initialName = initialFile.getName();
        int idx = initialName.lastIndexOf(".");
        if (idx > 0) {
            String ext = initialName.substring(idx);
            return fname + ext;
        }
        return fname;
    }


    private static Frame getParentFrame() {
        return IGV.hasInstance() ? IGV.getMainFrame() : null;
    }


}