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

package org.broad.igv.util;

import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.tools.IgvTools;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.Main;
import org.broad.tribble.util.ftp.FTPClient;
import org.junit.Ignore;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.util.concurrent.Future;

/**
 * @author jrobinso
 * @date Jul 28, 2010
 */
@Ignore
public class TestUtils {
    public static String DATA_DIR = "test/data";
    public static String dataFileName = DATA_DIR + "/genomes/hg18.unittest.genome";
    public static String AVAILABLE_FTP_URL = "ftp://ftp.broadinstitute.org/pub/igv/TEST/test.txt";
    public static String UNAVAILABLE_FTP_URL = "ftp://www.example.com/file.txt";
    //This is so ant can set the large data directory
    private static String LARGE_DATA_DIR_KEY = "LARGE_DATA_DIR";
    public static String LARGE_DATA_DIR = "test/largedata";

    static {
        LARGE_DATA_DIR = System.getProperty(LARGE_DATA_DIR_KEY, LARGE_DATA_DIR);
    }

    private static void setUpTestEnvironment() {
        Globals.setTesting(true);
        PreferenceManager.getInstance().setPrefsFile("testprefs.properties");
        Globals.READ_TIMEOUT = 60 * 1000;
        Globals.CONNECT_TIMEOUT = 60 * 1000;
        FTPClient.READ_TIMEOUT = 60 * 1000;
    }

    public static void setUpHeadless() {
        Globals.setHeadless(true);
        setUpTestEnvironment();
    }

    /**
     * This closes the IGV window.
     */
    public static void stopGUI() {
        try {
            IGV igv = IGV.getInstance();
        } catch (RuntimeException e) {
            return;
        }

        IGV.getMainFrame().setVisible(false);
        IGV.getMainFrame().dispose();
    }


    /**
     * Start GUI with default genome file
     *
     * @return
     * @throws IOException
     */
    public static IGV startGUI() throws IOException {
        return startGUI(dataFileName);
    }


    /**
     * Load a gui with the specified genome file.
     * No genome is loaded if null
     *
     * @param genomeFile
     * @return
     * @throws IOException
     */
    public static IGV startGUI(String genomeFile) throws IOException {
        boolean headless = checkHeadlessEnvironment();
        if (headless) {
            System.out.println("You are trying to start a GUI in a headless environment. Aborting test");
        }
        org.junit.Assume.assumeTrue(!headless);
        setUpTestEnvironment();
        IGV igv;
        //If IGV is already open, we get the instance.
        try {
            igv = IGV.getInstance();
            IGV.getMainFrame().setVisible(true);
            System.out.println("Using old IGV");
        } catch (RuntimeException e) {
            JFrame frame = new JFrame();
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            Main.open(frame);
            System.out.println("Started new IGV");
            igv = IGV.getInstance();
        }
        if (genomeFile != null) {
            igv.loadGenome(genomeFile, null);
        }
        return igv;
    }

    /**
     * Loads the session into IGV. This blocks until the session
     * is loaded.
     *
     * @param igv
     * @param sessionPath
     * @throws InterruptedException
     */
    public static void loadSession(IGV igv, String sessionPath) throws InterruptedException {
        Future worker = igv.doRestoreSession(sessionPath, null, false);

        while (!worker.isDone()) {
            Thread.sleep(100);
        }
    }

    public static boolean checkHeadlessEnvironment() {
        GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
        return ge.isHeadless();
    }

    /**
     * See TestUtils.createIndex(file, indexType, binSize)
     *
     * @param file
     * @throws IOException
     */
    public static void createIndex(String file) throws IOException {
        createIndex(file, IgvTools.LINEAR_INDEX, 1000);
    }

    /**
     * Destroys index file if it exists, and creates new one under
     * the specified parameters
     *
     * @param file
     * @param indexType
     * @param binSize
     * @throws IOException
     */
    public static void createIndex(String file, int indexType, int binSize) throws IOException {
        File indexFile = new File(file + ".idx");
        if (indexFile.exists()) {
            indexFile.delete();
        }
        (new IgvTools()).doIndex(file, indexType, binSize);
        indexFile.deleteOnExit();
    }

    /**
     * Load a test genome, do some test setup
     *
     * @return
     * @throws IOException
     */
    public static Genome loadGenome() throws IOException {
        final String genomeFile = dataFileName;
        return IgvTools.loadGenome(genomeFile, true);
    }

    public static void clearOutputDir() throws IOException {
        File outputDir = new File(DATA_DIR + "/out/");
        if (outputDir.isDirectory()) {
            File[] listFiles = outputDir.listFiles();
            for (File fi : listFiles) {
                fi.delete();
            }
        }
    }
}
