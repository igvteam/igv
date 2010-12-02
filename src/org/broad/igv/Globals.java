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

package org.broad.igv;

import org.apache.log4j.Logger;
import org.broad.igv.exceptions.DataLoadException;

import javax.swing.*;
import javax.swing.filechooser.FileSystemView;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;
import java.util.prefs.Preferences;
import java.util.regex.Pattern;

/**
 * User: jrobinso
 * Date: Feb 3, 2010
 */
public class Globals {

    private static Logger logger = Logger.getLogger(Globals.class);


    /**
     * CONSTANTS
     */
    final public static String CHR_ALL = "All";
    public static boolean headless = false;
    public static boolean suppressMessages = false;
    public static boolean batch = false;
    /**
     * Field description
     */
    final public static String SESSION_FILE_EXTENSION = ".xml";
    /**
     * GENOME ARCHIVE CONSTANTS
     */
    final public static String GENOME_FILE_EXTENSION = ".genome";
    final public static String ZIP_EXTENSION = ".zip";
    final public static String FASTA_GZIP_FILE_EXTENSION = ".gz";
    final public static String GENOME_ARCHIVE_PROPERTY_FILE_NAME = "property.txt";
    final public static String GENOME_ARCHIVE_ID_KEY = "id";
    final public static String GENOME_ARCHIVE_NAME_KEY = "name";
    final public static String GENOME_ARCHIVE_VERSION_KEY = "version";
    final public static String GENOME_ORDERED_KEY = "ordered";
    final public static String GENOME_GENETRACK_NAME = "geneTrackName";
    final public static String GENOME_ARCHIVE_CYTOBAND_FILE_KEY = "cytobandFile";
    final public static String GENOME_ARCHIVE_GENE_FILE_KEY = "geneFile";
    final public static String GENOME_ARCHIVE_SEQUENCE_FILE_LOCATION_KEY = "sequenceLocation";

    final static public Pattern commaPattern = Pattern.compile(",");
    final static public Pattern tabPattern = Pattern.compile("\t");
    final static public Pattern colonPattern = Pattern.compile(":");
    final static public Pattern dashPattern = Pattern.compile("-");
    final static public Pattern equalPattern = Pattern.compile("=");
    public static List emptyList = new ArrayList();
    public static String VERSION;
    public static String BUILD;
    public static String TIMESTAMP;
    public static final String NO_FEATURES_FOUND_WARNING = "No features were found in this file with chromosomes mapped to the current genome";
    public static double log2 = Math.log(2);

    /**
     * Field description
     */
    final public static boolean IS_WINDOWS =
            System.getProperty("os.name").toLowerCase().startsWith("windows");
    /**
     * Field description
     */
    final public static boolean IS_MAC =
            System.getProperty("os.name").toLowerCase().startsWith("mac");
    /**
     * Field description
     */
    final public static boolean IS_LINUX =
            System.getProperty("os.name").toLowerCase().startsWith("linux");
    final public static String GENOME_CACHE_FOLDER_NAME = "genomes";
    /**
     * Field description
     */
    private static File DEFAULT_USER_DIRECTORY;
    private static File DEFAULT_IGV_DIRECTORY;
    private static File GENOME_CACHE_DIRECTORY;
    private static File IGV_TEMP_DIRECTORY;
    private static File GENE_LIST_DIRECTORY;
    private static final String GENE_LIST_FOLDER_NAME = "lists";
    public static final String IGV_DIR_USERPREF = "igvDir";
    public static final String GENOME_CHR_ALIAS_FILE_KEY = "chrAliasFile";

    // Default user folder


    static {
        Properties properties = new Properties();
        try {
            properties.load(Globals.class.getResourceAsStream("/resources/about.properties"));
        }
        catch (IOException e) {
            logger.error("*** Error retrieving version and build information! ***", e);
        }
        VERSION = properties.getProperty("version", "???");
        BUILD = properties.getProperty("build", "???");
        TIMESTAMP = properties.getProperty("timestamp", "???");
    }

    public static void setHeadless(boolean bool) {
        headless = bool;
    }

    public static boolean isHeadless() {
        return headless;
    }

    public static void setSuppressMessages(boolean bool) {
        suppressMessages = bool;
    }

    public static boolean isSuppressMessages() {
        return suppressMessages;
    }

    public static String applicationString() {
        return "IGV Version " + VERSION + " (" + BUILD + ")" + TIMESTAMP;
    }

    public static String versionString() {
        return "<html>Version " + VERSION + " (" + BUILD + ")<br>" + TIMESTAMP;
    }

    public static synchronized File getUserDirectory() {
        if (DEFAULT_USER_DIRECTORY == null) {
            DEFAULT_USER_DIRECTORY = FileSystemView.getFileSystemView().getDefaultDirectory();
        }
        return DEFAULT_USER_DIRECTORY;
    }

    public static File getIgvDirectory() {

        // Hack for know Java bug
        if (System.getProperty("os.name").equals("Windows XP")) {
            try {
                Runtime.getRuntime().exec("attrib -r \"" + getUserDirectory().getAbsolutePath() + "\"");
            } catch (IOException e) {
                // Oh well, we tried

            }
        }

        if (DEFAULT_IGV_DIRECTORY == null)

        {

            // See if an override is stored in preferences.  Try to create a directory if it is.  If there is an
            // error (the log is likely not available yet) and try to use the standard directory
            try {
                Preferences prefs = Preferences.userNodeForPackage(Globals.class);
                String userDir = prefs.get(IGV_DIR_USERPREF, null);
                if (userDir != null) {
                    DEFAULT_IGV_DIRECTORY = new File(userDir);
                    if (!DEFAULT_IGV_DIRECTORY.exists()) {
                        DEFAULT_IGV_DIRECTORY = null;
                        prefs.remove(IGV_DIR_USERPREF);
                    }
                }
            } catch (Exception e) {
                Preferences prefs = Preferences.userNodeForPackage(Globals.class);
                prefs.remove(IGV_DIR_USERPREF);
                System.err.println("Error creating user directory");
                e.printStackTrace();
            }

            // No overide, try the default place
            if (DEFAULT_IGV_DIRECTORY == null) {
                String userHomeString = System.getProperty("user.home");
                File rootDir = new File(userHomeString);
                if (!(rootDir.exists() && canWrite(rootDir))) {
                    rootDir = getUserDirectory();
                }
                if (IS_MAC) {
                    DEFAULT_IGV_DIRECTORY = new File(rootDir, ".igv");
                } else {
                    DEFAULT_IGV_DIRECTORY = new File(rootDir, "igv");
                }
                if (!DEFAULT_IGV_DIRECTORY.exists()) {
                    try {
                        DEFAULT_IGV_DIRECTORY.mkdir();
                    }
                    catch (Exception e) {
                        System.err.println("Error creating user directory");
                        e.printStackTrace();
                    }
                }
            }

            // The IGV directory either doesn't exist or isn't writeable.  This situation can arise with Windows Vista
            // and Windows 7 due to a Java bug (http://bugs.sun.com/view_bug.do?bug_id=4787931)
            if (!(DEFAULT_IGV_DIRECTORY.exists() && DEFAULT_IGV_DIRECTORY.canRead() && canWrite(DEFAULT_IGV_DIRECTORY))) {
                if (isHeadless() || isSuppressMessages()) {
                    System.err.println("Cannot write to igv directory: " + DEFAULT_IGV_DIRECTORY.getAbsolutePath());
                    DEFAULT_IGV_DIRECTORY = (new File(".")).getParentFile();
                } else {
                    int option = JOptionPane.showConfirmDialog(null,
                            "<html>The default IGV directory (" + DEFAULT_IGV_DIRECTORY + ") " +
                                    "cannot be accessed.  Click Yes to choose a new folder or No to exit.<br>" +
                                    "This folder will be used to store user preferences and cached genomes.",
                            "IGV Directory Error", JOptionPane.YES_NO_OPTION);

                    if (option == JOptionPane.YES_OPTION) {
                        final JFileChooser fc = new JFileChooser();
                        fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
                        int retValue = fc.showOpenDialog(null);
                        if (retValue == JFileChooser.APPROVE_OPTION) {
                            DEFAULT_IGV_DIRECTORY = fc.getSelectedFile();
                            Preferences prefs = Preferences.userNodeForPackage(Globals.class);
                            prefs.put(IGV_DIR_USERPREF, DEFAULT_IGV_DIRECTORY.getAbsolutePath());
                        }
                    }
                }
            }


            if (!DEFAULT_IGV_DIRECTORY.canRead()) {
                throw new DataLoadException("Cannot read from user directory", DEFAULT_IGV_DIRECTORY.getAbsolutePath());
            } else if (!canWrite(DEFAULT_IGV_DIRECTORY)) {
                throw new DataLoadException("Cannot write to user directory", DEFAULT_IGV_DIRECTORY.getAbsolutePath());
            }
        }


        return DEFAULT_IGV_DIRECTORY;
    }

    private static boolean canWrite(File directory) {
        // There are bugs in Java window (targe fix is Java 7).  The only way to know for sure is to try to write something
        if (IS_WINDOWS) {
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


    public static File getGenomeCacheDirectory() {
        if (GENOME_CACHE_DIRECTORY == null) {

            //Create the Genome Cache
            GENOME_CACHE_DIRECTORY = new File(getIgvDirectory(), GENOME_CACHE_FOLDER_NAME);
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

            //Create the Genome Cache
            GENE_LIST_DIRECTORY = new File(getIgvDirectory(), GENE_LIST_FOLDER_NAME);
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

    public static File getTempDirectory() {

        if (IGV_TEMP_DIRECTORY == null) {
            String path = System.getProperty("java.io.tmpdir");
            IGV_TEMP_DIRECTORY = new File(path.toString());
            if (!IGV_TEMP_DIRECTORY.exists()) {
                IGV_TEMP_DIRECTORY.mkdir();
            }
            if (!IGV_TEMP_DIRECTORY.canRead()) {
                throw new DataLoadException("Cannot read from user directory", IGV_TEMP_DIRECTORY.getAbsolutePath());
            } else if (!IGV_TEMP_DIRECTORY.canWrite()) {
                throw new DataLoadException("Cannot write to user directory", IGV_TEMP_DIRECTORY.getAbsolutePath());
            }
        }
        return IGV_TEMP_DIRECTORY;
    }
}
