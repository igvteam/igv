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

package org.broad.igv.ui.util;

import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.ui.IGV;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.io.FilenameFilter;
import java.lang.reflect.Method;

/**
 * @author jrobinso
 * @date Nov 17, 2010
 */
public class FileDialogUtils {

    public static int LOAD = FileDialog.LOAD;
    public static int SAVE = FileDialog.SAVE;

    private static Logger log = Logger.getLogger(FileDialogUtils.class);


    public static File chooseFile(String title, File initialDirectory, int mode) {
        return chooseFile(title, initialDirectory, null, mode);
    }

    public static File chooseFile(String title) {
        return chooseFile(title, DirectoryManager.getUserDirectory(), null, FileDialog.LOAD);
    }

    public static File chooseFile(String title, File initialDirectory, File initialFile, int mode) {
        return chooseFile(title, initialDirectory, initialFile, null, JFileChooser.FILES_ONLY, mode);
    }

    private static File chooseFile(String title, File initialDirectory, File initialFile, FilenameFilter filter,
                                   int directoriesMode, int mode) {

        if (initialDirectory == null && initialFile != null) {
            initialDirectory = initialFile.getParentFile();
        }
        // Strip off parent directory
        if (initialFile != null) initialFile = new File(initialFile.getName());
        return chooseNative(title, initialDirectory, initialFile, filter, directoriesMode, mode);

    }

    public static File chooseDirectory(String title, File initialDirectory) {
        return chooseNative(title, initialDirectory, null, null, JFileChooser.DIRECTORIES_ONLY, LOAD);
    }

    public static File[] chooseMultiple(String title, File initialDirectory, final FilenameFilter filter) {

        File[] files = null;
        FileDialog fd = getNativeChooser(title, initialDirectory, null, filter, JFileChooser.FILES_ONLY, LOAD);
        if (fd.isMultipleMode()) {
            fd.setVisible(true);
            files = fd.getFiles();
        }
        return files;
    }

    private static FileDialog getNativeChooser(String title, File initialDirectory, File initialFile, FilenameFilter filter, int directoryMode, int mode) {
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

        if (mode == LOAD && !directories) {
            setMultipleMode(fd, true);
        }
        return fd;
    }


    private static File chooseNative(String title, File initialDirectory, File initialFile, FilenameFilter filter,
                                     int directoryMode, int mode) {

        FileDialog fd = getNativeChooser(title, initialDirectory, initialFile, filter, directoryMode, mode);
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

    /**
     * Reflectively call FileDialog.setMultipleMode.
     * Does nothing if method not available
     *
     * @param fd
     * @param b
     * @return true if call was successful, false if not
     */
    private static boolean setMultipleMode(FileDialog fd, boolean b) {
        try {
            Method method = FileDialog.class.getMethod("setMultipleMode", boolean.class);
            method.invoke(fd, b);
            return true;
        } catch (Exception e) {
            return false;
        }
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