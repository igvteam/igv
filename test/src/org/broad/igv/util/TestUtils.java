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
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.tools.IgvTools;
import org.broad.tribble.util.ftp.FTPClient;

import java.io.File;
import java.io.IOException;

/**
 * @author jrobinso
 * @date Jul 28, 2010
 */
public class TestUtils {
    public static String DATA_DIR = "test/data";
    static String dataFileName = DATA_DIR + "/genomes/hg18.unittest.genome";
    public static String AVAILABLE_FTP_URL = "ftp://ftp.broadinstitute.org/pub/igv/TEST/test.txt";
    public static String UNAVAILABLE_FTP_URL = "ftp://www.example.com/file.txt";
    public static String LARGE_DATA_DIR = "test/largedata";


//    public static Genome loadGenome(String genomeId) throws IOException {
//        Globals.setHeadless(true);
//        return IgvTools.loadGenome(genomeId, true);
//
//    }

    public static void setUpTestEnv() {
        Globals.setHeadless(true);
        Globals.READ_TIMEOUT = 30 * 1000;
        Globals.CONNECT_TIMEOUT = 30 * 1000;
        FTPClient.READ_TIMEOUT = 30 * 1000;
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
