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

import org.apache.log4j.Logger;
import org.broad.igv.renderer.SequenceRenderer;

import java.io.IOException;
import java.util.*;
import java.util.List;
import java.util.regex.Pattern;

/**
 * User: jrobinso
 * Date: Feb 3, 2010
 */
public class Globals {

    public static final int DESIGN_DPI = 96;
    private static Logger log = Logger.getLogger(Globals.class);


    /**
     * CONSTANTS
     */
    public static final String CHR_ALL = "All";
    public static final String TRACK_NAME_ATTRIBUTE = "NAME";
    public static final String TRACK_DATA_FILE_ATTRIBUTE = "DATA FILE";
    public static final String TRACK_DATA_TYPE_ATTRIBUTE = "DATA TYPE";

    private static boolean headless = false;
    private static boolean suppressMessages = false;
    private static boolean batch = false;
    private static boolean testing = false;
    public static int CONNECT_TIMEOUT = 20000;        // 20 seconds
    public static int READ_TIMEOUT = 1000 * 3 * 60;   // 3 minutes

    /**
     * Field description
     */
    final public static String SESSION_FILE_EXTENSION = ".xml";
    /**
     * GENOME ARCHIVE CONSTANTS
     */
    final public static String GENOME_FILE_EXTENSION = ".genome";
    final public static String ZIP_EXTENSION = ".zip";
    final public static String GZIP_FILE_EXTENSION = ".gz";
    final public static String GENOME_ARCHIVE_PROPERTY_FILE_NAME = "property.txt";
    final public static String GENOME_ARCHIVE_ID_KEY = "id";
    final public static String GENOME_ARCHIVE_NAME_KEY = "name";
    final public static String GENOME_ARCHIVE_VERSION_KEY = "version";
    final public static String GENOME_ORDERED_KEY = "ordered";
    final public static String GENOME_GENETRACK_NAME = "geneTrackName";
    final public static String GENOME_URL_KEY = "url";
    final public static String GENOME_ARCHIVE_CYTOBAND_FILE_KEY = "cytobandFile";
    final public static String GENOME_ARCHIVE_GENE_FILE_KEY = "geneFile";
    final public static String GENOME_ARCHIVE_SEQUENCE_FILE_LOCATION_KEY = "sequenceLocation";

    /**
     * Whether the sequenceLocation has been modified from the version of the .genome
     * file on the server
     */
    public static final String GENOME_ARCHIVE_CUSTOM_SEQUENCE_LOCATION_KEY = "customSequenceLocation";
    public static final String GENOME_CHR_ALIAS_FILE_KEY = "chrAliasFile";

    // Default user folder

    final static public Pattern commaPattern = Pattern.compile(",");
    final static public Pattern tabPattern = Pattern.compile("\t");
    final static public Pattern colonPattern = Pattern.compile(":");
    final static public Pattern dashPattern = Pattern.compile("-");
    final static public Pattern equalPattern = Pattern.compile("=");
    final static public Pattern semicolonPattern = Pattern.compile(";");
    final static public Pattern singleTabMultiSpacePattern = Pattern.compile("\t|( +)");
    final static public Pattern forwardSlashPattern = Pattern.compile("/");
    final static public Pattern whitespacePattern = Pattern.compile("\\s+");


    public static List emptyList = new ArrayList();
    public static String VERSION;
    public static String BUILD;
    public static String TIMESTAMP;
    public static double log2 = Math.log(2);


    final public static boolean IS_WINDOWS =
            System.getProperty("os.name").toLowerCase().startsWith("windows");
    final public static boolean IS_MAC =
            System.getProperty("os.name").toLowerCase().startsWith("mac");

    final public static boolean IS_LINUX =
            System.getProperty("os.name").toLowerCase().startsWith("linux");

    final public static boolean IS_JWS =
            System.getProperty("webstart.version", null) != null || System.getProperty("javawebstart.version", null) != null;

    public static final String JAVA_VERSION_STRING = "java.version";

    //Location of bedtools executable
    //Note: It is recommended you use an absolute path here.
    //System paths can be finnicky and vary depending on how IGV is launched
    //However, the path environment variable will be checked if the executable
    //is named rather than the full path given
    public static String BEDtoolsPath = "/usr/local/bin/bedtools"; //"bedtools"
    public static boolean toolsMenuEnabled = false;
    public static boolean development;

    public static String versionURL = "http://www.broadinstitute.org/igv/projects/current/version.txt";
    public static String downloadURL = "http://www.broadinstitute.org/igv/download";
    static {
        Properties properties = new Properties();
        try {
            properties.load(Globals.class.getResourceAsStream("/resources/about.properties"));
        } catch (IOException e) {
            log.error("*** Error retrieving version and build information! ***", e);
        }
        VERSION = properties.getProperty("version", "???");
        BUILD = properties.getProperty("build", "???");
        TIMESTAMP = properties.getProperty("timestamp", "???");
        BEDtoolsPath = System.getProperty("BEDtoolsPath", BEDtoolsPath);

         //Runtime property overrides compile-time property, if both exist.
        //If neither exist we default to false
        final String prodProperty = System.getProperty("development", properties.getProperty("development", "false"));
        development = Boolean.parseBoolean(prodProperty);
        if(development){
            log.warn("Development mode is enabled");
        }

    }

    public static void setHeadless(boolean bool) {
        headless = bool;
    }

    public static boolean isHeadless() {
        return headless;
    }

    public static void setTesting(boolean testing) {
        Globals.testing = testing;
    }

    public static boolean isTesting() {
        return testing;
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

    public static boolean isDevelopment() {
        return development;
    }


    public static boolean isBatch() {
        return batch;
    }

    public static void setBatch(boolean batch) {
        Globals.batch = batch;
    }


    /**
     * Checks whether the current JVM is the minimum specified version
     * or higher. Only compares up to as many characters as
     * in {@code minVersion}
     *
     * @param minVersion
     * @return
     */
    public static boolean checkJavaVersion(String minVersion) {
        String curVersion = System.getProperty(JAVA_VERSION_STRING);
        if (curVersion.length() >= minVersion.length()) {
            curVersion = curVersion.substring(0, minVersion.length());
        }
        return curVersion.compareTo(minVersion) >= 0;
    }

    /**
     * Return the URL to fetch the current IGV version (note:  not Java version)
     * @return
     */
    public static String getVersionURL() {
        return versionURL;
    }
}
