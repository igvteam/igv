/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv;

import org.apache.log4j.Logger;

import java.awt.*;
import java.io.IOException;
import java.util.*;
import java.util.List;
import java.util.regex.Pattern;

/**
 * User: jrobinso
 * Date: Feb 3, 2010
 */
public class Globals {

    private static Logger log = Logger.getLogger(Globals.class);


    /**
     * CONSTANTS
     */
    final public static String CHR_ALL = "All";
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
    public static final String GENOME_CHR_ALIAS_FILE_KEY = "chrAliasFile";
    public static final String DEFAULT_GENOME = "hg18";

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

    public static final String JAVA_VERSION_STRING = "java.version";
    public static Map<Character, Color> nucleotideColors;

    //Location of bedtools executable
    //Note: It is recommended you use an absolute path here.
    //System paths can be finnicky and vary depending on how IGV is launched
    //However, the path environment variable will be checked if the executable
    //is named rather than the full path given
    public static String BEDtoolsPath = "/usr/local/bin/bedtools"; //"bedtools"
    public static boolean toolsMenuEnabled = false;

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

        nucleotideColors = new HashMap();
        nucleotideColors.put('A', Color.GREEN);
        nucleotideColors.put('a', Color.GREEN);
        nucleotideColors.put('C', Color.BLUE);
        nucleotideColors.put('c', Color.BLUE);
        nucleotideColors.put('T', Color.RED);
        nucleotideColors.put('t', Color.RED);
        nucleotideColors.put('G', new Color(209, 113, 5));
        nucleotideColors.put('g', new Color(209, 113, 5));
        nucleotideColors.put('N', Color.gray.brighter());
        nucleotideColors.put('n', Color.gray.brighter());

        BEDtoolsPath = System.getProperty("BEDtoolsPath", BEDtoolsPath);

        toolsMenuEnabled = Boolean.parseBoolean(System.getProperty("enable.tools", "false"));
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


    public static boolean isBatch() {
        return batch;
    }

    public static void setBatch(boolean batch) {
        Globals.batch = batch;
    }

}
