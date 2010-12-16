/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

import org.broad.igv.Globals;
import org.broad.igv.ui.IGVMainFrame;

import javax.swing.*;
import java.awt.*;
import java.io.File;

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
        return chooseFile(title, Globals.getUserDirectory(), null, FileDialog.LOAD);
    }

    public static File chooseFile(String title, File initialDirectory, File initialFile, int mode) {

        if(initialDirectory == null && initialFile != null) {
            initialDirectory = initialFile.getParentFile();
        }

        // Strip off parent directory
        if(initialFile != null) initialFile = new File(initialFile.getName());

        if (Globals.IS_MAC) {
            return chooseNative(title, initialDirectory, initialFile, false, mode);
        } else {
            return chooseSwing(title, initialDirectory, initialFile, false, mode);
        }
    }

    public static File chooseDirectory(String title, File initialDirectory) {

        if (Globals.IS_MAC) {
            return chooseNative(title, initialDirectory, null, true, LOAD);
        } else {
            return chooseSwing(title, initialDirectory, null, true, LOAD);
        }
    }

    private static File chooseNative(String title, File initialDirectory, File initialFile, boolean directories, int mode) {

        System.setProperty("apple.awt.fileDialogForDirectories", String.valueOf(directories));
        FileDialog fd = new FileDialog(IGVMainFrame.getInstance(), title);
        if (initialDirectory != null) {
            fd.setDirectory(initialDirectory.getAbsolutePath());
        }
        if (initialFile != null) {
            fd.setFile(initialFile.getName());
        }
        fd.setModal(true);
        fd.setMode(mode);
        fd.setVisible(true);

        String file = fd.getFile();
        String directory = fd.getDirectory();
        if (file != null && directory != null) {
            return new File(directory, file);
        } else {
            return null;
        }
    }


    private static File chooseSwing(String title, File initialDirectory, File initialFile, boolean directories, int mode) {


        JFileChooser fileChooser = new JFileChooser();
        if (initialDirectory != null) {
            fileChooser.setCurrentDirectory(initialDirectory);
        }
         if (initialFile != null) {
            fileChooser.setSelectedFile(initialFile);
        }
        fileChooser.setDialogTitle(title);
        fileChooser.setFileSelectionMode(directories ? JFileChooser.DIRECTORIES_ONLY : JFileChooser.FILES_ONLY);

        boolean approve = false;
        if (mode == LOAD) {
            approve = fileChooser.showOpenDialog(IGVMainFrame.getInstance()) == JFileChooser.APPROVE_OPTION;
        } else {
            approve = fileChooser.showSaveDialog(IGVMainFrame.getInstance()) == JFileChooser.APPROVE_OPTION;
        }

        if (approve) {
            return fileChooser.getSelectedFile();
        } else {
            return null;
        }

    }

}