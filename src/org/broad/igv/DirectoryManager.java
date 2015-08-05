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

package org.broad.igv;

import org.apache.commons.io.FileUtils;
import org.apache.log4j.*;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.RuntimeUtils;

import javax.swing.*;
import javax.swing.filechooser.FileSystemView;
import java.io.File;
import java.io.IOException;
import java.util.prefs.Preferences;

/**
 * @author Jim Robinson
 * @date 3/19/12
 */
public class DirectoryManager {

    private static Logger log = Logger.getLogger(DirectoryManager.class);

    private static File USER_HOME;
    private static File USER_DIRECTORY;    // FileSystemView.getFileSystemView().getDefaultDirectory();
    private static File IGV_DIRECTORY;     // The IGV application directory
    private static File GENOME_CACHE_DIRECTORY;
    private static File GENE_LIST_DIRECTORY;
    private static File BAM_CACHE_DIRECTORY;
    final public static String IGV_DIR_USERPREF = "igvDir";


    private static File getUserHome() {
        if (USER_HOME == null) {
            String userHomeString = System.getProperty("user.home");
            USER_HOME = new File(userHomeString);
        }
        return USER_HOME;
    }

    /**
     * The user directory.  On Mac and Linux this should be the user home directory.  On Windows platforms this
     * is the "My Documents" directory.
     */
    public static synchronized File getUserDirectory() {
        if (USER_DIRECTORY == null) {
            log.info("Fetching user directory... ");
            USER_DIRECTORY = FileSystemView.getFileSystemView().getDefaultDirectory();
            //Mostly for testing, in some environments USER_DIRECTORY can be null
            if (USER_DIRECTORY == null) {
                USER_DIRECTORY = getUserHome();
            }
        }
        return USER_DIRECTORY;
    }


    public static File getIgvDirectory() {

        if (IGV_DIRECTORY == null) {

            // Hack for known Java / Windows bug.   Attempt to remvoe (possible) read-only bit from user directory
            if (System.getProperty("os.name").startsWith("Windows")) {
                try {
                    Runtime.getRuntime().exec("attrib -r \"" + getUserDirectory().getAbsolutePath() + "\"");
                } catch (Exception e) {
                    // We tried
                }
            }

            IGV_DIRECTORY = getIgvDirectoryOverride();

            // If still null, try the default place
            if (IGV_DIRECTORY == null) {
                File rootDir = getUserHome();
                IGV_DIRECTORY = new File(rootDir, "igv");

                if (!IGV_DIRECTORY.exists()) {
                    // See if a pre-2.1 release directory exists, if so copy it
                    File legacyDirectory = null;
                    try {
                        legacyDirectory = getLegacyIGVDirectory();
                        if (legacyDirectory.exists()) {
                            log.info("Copying " + legacyDirectory + " => " + IGV_DIRECTORY);
                            FileUtils.copyDirectory(legacyDirectory, IGV_DIRECTORY);
                        }
                    } catch (IOException e) {
                        log.error("Error copying igv directory " + legacyDirectory + " => " + IGV_DIRECTORY, e);
                    }
                }

                if (!IGV_DIRECTORY.exists()) {
                    try {
                        boolean wasSuccessful = IGV_DIRECTORY.mkdir();
                        if (!wasSuccessful) {
                            log.error("Failed to create user directory!");
                            IGV_DIRECTORY = null;
                        }
                    } catch (Exception e) {
                        log.error("Error creating igv directory", e);
                    }
                }
            }


            // The IGV directory either doesn't exist or isn't writeable.  This situation can arise with Windows Vista
            // and Windows 7 due to a Java bug (http://bugs.sun.com/view_bug.do?bug_id=4787931)
            if (IGV_DIRECTORY == null || !IGV_DIRECTORY.exists() || !canWrite(IGV_DIRECTORY)) {
                if (Globals.isHeadless() || Globals.isSuppressMessages()) {
                    System.err.println("Cannot write to igv directory: " + IGV_DIRECTORY.getAbsolutePath());
                    IGV_DIRECTORY = (new File(".")).getParentFile();
                } else {
                    int option = JOptionPane.showConfirmDialog(null,
                            "<html>The default IGV directory (" + IGV_DIRECTORY + ") " +
                                    "cannot be accessed.  Click Yes to choose a new folder or No to exit.<br>" +
                                    "This folder will be used to create the 'igv' directory",
                            "IGV Directory Error", JOptionPane.YES_NO_OPTION);

                    if (option == JOptionPane.YES_OPTION) {
                        File parentDirectory = FileDialogUtils.chooseDirectory("Select a location for the igv directory", null);
                        if (parentDirectory != null) {
                            IGV_DIRECTORY = new File(parentDirectory, "igv");
                            IGV_DIRECTORY.mkdir();
                            Preferences prefs = Preferences.userNodeForPackage(Globals.class);
                            prefs.put(IGV_DIR_USERPREF, IGV_DIRECTORY.getAbsolutePath());
                        }
                    }
                }
            }


            if (IGV_DIRECTORY == null || !IGV_DIRECTORY.canRead()) {
                throw new DataLoadException("Cannot read from user directory", IGV_DIRECTORY.getAbsolutePath());
            } else if (!canWrite(IGV_DIRECTORY)) {
                throw new DataLoadException("Cannot write to user directory", IGV_DIRECTORY.getAbsolutePath());
            }

            log.info("IGV Directory: " + IGV_DIRECTORY.getAbsolutePath());
        }
        return IGV_DIRECTORY;
    }

    private static File getIgvDirectoryOverride() {
        Preferences userPrefs = null;
        File override = null;
        try {
            // See if an override location has been specified.  This is stored with the Java Preferences API
            userPrefs = Preferences.userNodeForPackage(Globals.class);
            String userDir = userPrefs.get(IGV_DIR_USERPREF, null);
            if (userDir != null) {
                override = new File(userDir);
                if (!override.exists()) {
                    override = null;
                    userPrefs.remove(IGV_DIR_USERPREF);
                }
            }
        } catch (Exception e) {
            userPrefs.remove(IGV_DIR_USERPREF);
            override = null;
            log.error("Error creating user directory", e);
        }
        return override;
    }


    public static File getGenomeCacheDirectory() {
        if (GENOME_CACHE_DIRECTORY == null) {

            //Create the Genome Cache
            GENOME_CACHE_DIRECTORY = new File(getIgvDirectory(), "genomes");
            if (!GENOME_CACHE_DIRECTORY.exists()) {
                GENOME_CACHE_DIRECTORY.mkdir();
            }
            if (!GENOME_CACHE_DIRECTORY.canRead()) {
                throw new DataLoadException("Cannot read from user directory", GENOME_CACHE_DIRECTORY.getAbsolutePath());
            } else if (!GENOME_CACHE_DIRECTORY.canWrite()) {
                throw new DataLoadException("Cannot write to user directory", GENOME_CACHE_DIRECTORY.getAbsolutePath());
            }
        }
        return GENOME_CACHE_DIRECTORY;
    }

    public static File getGeneListDirectory() {
        if (GENE_LIST_DIRECTORY == null) {
            GENE_LIST_DIRECTORY = new File(getIgvDirectory(), "lists");
            if (!GENE_LIST_DIRECTORY.exists()) {
                GENE_LIST_DIRECTORY.mkdir();
            }
            if (!GENE_LIST_DIRECTORY.canRead()) {
                throw new DataLoadException("Cannot read from user directory", GENE_LIST_DIRECTORY.getAbsolutePath());
            } else if (!GENE_LIST_DIRECTORY.canWrite()) {
                throw new DataLoadException("Cannot write to user directory", GENE_LIST_DIRECTORY.getAbsolutePath());
            }
        }
        return GENE_LIST_DIRECTORY;
    }

    public static synchronized File getCacheDirectory() {
        if (BAM_CACHE_DIRECTORY == null) {
            File defaultDir = getIgvDirectory();
            if (defaultDir.exists()) {
                BAM_CACHE_DIRECTORY = new File(defaultDir, "bam");
                if (!BAM_CACHE_DIRECTORY.exists()) {
                    BAM_CACHE_DIRECTORY.mkdir();
                }
            }
        }
        return BAM_CACHE_DIRECTORY;
    }

    public static synchronized File getSamDirectory() {

        File samDir = new File(DirectoryManager.getIgvDirectory(), "sam");
        if (!samDir.exists()) {
            samDir.mkdir();
        }
        return samDir;

    }

    public static synchronized File getLogFile() throws IOException {

        File logFile = new File(getIgvDirectory(), "igv.log");
        if (!logFile.exists()) {
            logFile.createNewFile();
        }
        return logFile;

    }


    /**
     * Return the user preferences property file  ("~/igv/prefs.properties").
     *
     * @return
     */
    public static synchronized File getPreferencesFile() {

        File igvDirectoy = getIgvDirectory();
        File igvPropertyFile = new File(igvDirectoy, "prefs.properties");

        // If the property file doesn't exist, try the "legacy" location.  This should only make a difference for Macs
        if (!igvPropertyFile.exists()) {
            File oldFile = getLegacyPreferencesFile();
            if (oldFile.exists()) {
                try {
                    // Copy contents to new location.  Leave oldFile in place for compatibility with earlier IGV releases
                    FileUtils.copyFile(oldFile, igvPropertyFile);
                } catch (IOException e) {
                    log.error("Error copy property file from: " + oldFile + " to: " + igvPropertyFile, e);
                }
            }
        }

        if (!igvPropertyFile.exists()) {
            try {
                igvPropertyFile.createNewFile();
            } catch (IOException e) {
                log.error("Could not create property file: " + igvPropertyFile, e);
            }
        }

        return igvPropertyFile;
    }


    /**
     * Move the "igv" directory to a new location, copying all contents.  Returns True if the directory
     * is successfully moved, irrespective of any errors that might occur later (e.g. when attempting to
     * remove the old directory).
     *
     * @param newIGVDirectory
     * @return True if the directory is successfully moved, false otherwise
     */

    public static Boolean moveIGVDirectory(final File newIGVDirectory) {

        if (newIGVDirectory.equals(IGV_DIRECTORY)) {
            return false; // Nothing to do
        }

        if (IGV_DIRECTORY != null && IGV_DIRECTORY.exists()) {

            File oldDirectory = IGV_DIRECTORY;

            try {
                log.info("Moving igv directory from " + oldDirectory.getParent() + " to " +
                        newIGVDirectory.getAbsolutePath());
                FileUtils.copyDirectory(IGV_DIRECTORY, newIGVDirectory);
                IGV_DIRECTORY = newIGVDirectory;

                // Store location of new directory in Java preferences node (not pref.properties)
                Preferences prefs = Preferences.userNodeForPackage(Globals.class);
                prefs.put(IGV_DIR_USERPREF, newIGVDirectory.getAbsolutePath());

                // Update preference manager with new file location
                PreferenceManager.getInstance().setPrefsFile(getPreferencesFile().getAbsolutePath());

            } catch (IOException e) {
                log.error("Error copying IGV directory", e);
                MessageUtils.showMessage("<html>Error moving IGV directory:<br/>&nbsp;nbsp;" + e.getMessage());
                return false;
            }

            // Restart the log
            LogManager.shutdown();
            initializeLog();

            // Try to delete the old directory
            try {
                deleteDirectory(oldDirectory);
            } catch (IOException e) {
                log.error("An error was encountered deleting the previous IGV directory", e);
                MessageUtils.showMessage("<html>An error was encountered deleting the previous IGV directory (" +
                        e.getMessage() + "):<br>&nbsp;nbsp;nbsp;" + oldDirectory.getAbsolutePath() +
                        "<br>Remaining files should be manually deleted.");
            }

        }

        GENOME_CACHE_DIRECTORY = null;
        GENE_LIST_DIRECTORY = null;
        BAM_CACHE_DIRECTORY = null;
        return true;

    }

    /**
     * Delete the directory and all contents recursively.  The apache FileUtils is hanging on Linux.
     *
     * @param oldDirectory
     * @throws IOException
     */
    private static void deleteDirectory(File oldDirectory) throws IOException {
        if (Globals.IS_LINUX || Globals.IS_MAC) {
            //System.out.println("Deleting: " + oldDirectory);
            String[] cmd = new String[]{"rm", "-rf", oldDirectory.getAbsolutePath()};
            String result = RuntimeUtils.executeShellCommand(cmd, null, null);
            if (result != null && result.trim().length() > 0) {
                log.info("Response from 'rm -rf': " + result);
            }
        } else {
            // The apache commons FileUtils is not working reliably
            org.broad.igv.util.FileUtils.deleteDir(oldDirectory);
        }
    }


    /**
     * Return the pre 2.1 release properties file.  This may or may not exist.
     *
     * @return
     */
    private static synchronized File getLegacyPreferencesFile() {
        File rootDir = getLegacyIGVDirectory();
        return new File(rootDir, "prefs.properties");
    }

    /**
     * REturn the pre 2.1 release IGV directory.  This differs from the current location only for Macs.
     *
     * @return
     */
    private static File getLegacyIGVDirectory() {
        File rootDir = getUserHome();
        if (Globals.IS_MAC) {
            rootDir = new File(rootDir, ".igv");
        } else {
            rootDir = new File(rootDir, "igv");
        }
        return rootDir;
    }

    private static boolean canWrite(File directory) {
        // There are bugs in the Windows Java JVM that can cause user directories to be non-writable (target fix is
        // Java 7).  The only way to know if the directory is writable for sure is to try to write something.
        if (Globals.IS_WINDOWS) {
            File testFile = null;
            try {
                testFile = new File(directory, "igv332415dsfjdsklt.testfile");
                if (testFile.exists()) {
                    testFile.delete();
                }
                testFile.deleteOnExit();
                testFile.createNewFile();
                return testFile.exists();
            } catch (IOException e) {
                return false;
            } finally {
                if (testFile.exists()) {
                    testFile.delete();
                }
            }
        } else {
            return directory.canWrite();
        }

    }

    public static void initializeLog() {

        Logger logger = Logger.getRootLogger();
        PatternLayout layout = new PatternLayout();
        layout.setConversionPattern("%p [%d{ISO8601}] [%F:%L]  %m%n");

        //Temp appender just used while initializing log
        //If something goes wrong here we want a log of it
        ConsoleAppender tempAppender = new ConsoleAppender(layout, "System.err");
        logger.addAppender(tempAppender);

        // Create a log file that is ready to have text appended to it
        try {
            File logFile = getLogFile();
            RollingFileAppender appender = new RollingFileAppender();
            appender.setName("IGV_ROLLING_APPENDER");
            appender.setFile(logFile.getAbsolutePath());
            //Will actually inherit from root logger
            appender.setThreshold(Level.ALL);
            appender.setMaxFileSize("1000KB");
            appender.setMaxBackupIndex(1);
            appender.setLayout(layout);
            appender.setAppend(true);
            appender.activateOptions();
            logger.addAppender(appender);
            logger.removeAppender(tempAppender);
        } catch (IOException e) {
            // Can't create log file, just log to console as set in log4j.properties
            log.error("Error creating log file", e);
        }
    }
}
