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

package org.broad.igv.tools;

import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.tools.parsers.DataConsumer;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.LineReader;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.assertEquals;


public class CoverageCounterTest {

    static PreferenceManager preferenceManager;
    static boolean useByteRange;
    static Genome genome;

    @BeforeClass
    public static void setUpClass() throws Exception {
        genome = TestUtils.loadGenome();
        preferenceManager = PreferenceManager.getInstance();
        useByteRange = preferenceManager.getAsBoolean(PreferenceManager.USE_BYTE_RANGE);
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
        preferenceManager.put(PreferenceManager.USE_BYTE_RANGE, useByteRange);
        //TestUtils.clearOutputDir();
    }

    @Before
    public void setUp() {
        preferenceManager.put(PreferenceManager.USE_BYTE_RANGE, true);
    }

    /**
     * Test the "mapping quality" flag.  Also indirectly tests the query parameters.
     */
    @Test
    public void testMappingQualityFlag() throws IOException {
        String bamURL = "http://www.broadinstitute.org/igvdata/1KG/pilot2Bams/NA12878.SLX.bam";
        String options = "m=30,q@1:16731624-16731624";
        File wigFile = null;
        Genome genome = null;
        int windowSize = 1;

        TestDataConsumer dc = new TestDataConsumer();

        CoverageCounter cc = new CoverageCounter(bamURL, dc, windowSize, 0, wigFile, genome, options);

        cc.parse();

        String totalCount = dc.attributes.get("totalCount");

        assertEquals("19", totalCount);

    }

    /*
    Test whether we count the number of features correctly.
     */
    @Test
    public void testCountStrand() throws Exception {
        String ifile = TestUtils.DATA_DIR + "/bed/Unigene.sample.bed";


        int expTot = 71;
        int windowSize = 1;

        File wigFile = new File(TestUtils.DATA_DIR + "/out", "testCountStrand.wig");
        Genome genome = this.genome;


        String[] options = new String[]{"sc=0x01", "sc=0x04", "sc=0x08", "sc=0xC", "sc=0x1C", "sc=0x05", "sc=0x02"};
        int[] expectedTotal = new int[]{expTot, expTot, expTot, expTot, expTot, expTot, 0};

        for (int ii = 0; ii < options.length; ii++) {
            TestDataConsumer dc = new TestDataConsumer();
            CoverageCounter cc = new CoverageCounter(ifile, dc, windowSize, 0, wigFile, genome, options[ii]);
            cc.parse();

            String totalCount = dc.attributes.get("totalCount");
            int totCount = Integer.valueOf(totalCount);
            assertEquals(expectedTotal[ii], totCount);
        }
    }

    /*
    Simple test, output all 3 data columns and check that they add up properly (both = positive + negative)
    */
    @Test
    public void testStrandsConsistent() throws Exception {
        String ifile = TestUtils.DATA_DIR + "/bed/Unigene.sample.bed";
        int[] windowSizes = new int[]{1, 10, 50, 101, 500, 999};
        int expected_cols = 4;

        File wigFile = new File(TestUtils.DATA_DIR + "/out", "testStrandsConsistent.wig");
        Genome genome = this.genome;


        for (int ii = 0; ii < windowSizes.length; ii++) {
            TestDataConsumer dc = new TestDataConsumer();
            CoverageCounter cc = new CoverageCounter(ifile, dc, windowSizes[ii], 0, wigFile, genome, "sc=" + (1 + 4 + 8));
            cc.parse();

            InputStream is = new FileInputStream(wigFile);
            LineReader data = new AsciiLineReader(is);
            //Skip header
            String line = data.readLine();
            while ((line = data.readLine()) != null) {
                String[] tokens = line.split("\\s+");
                assertEquals(expected_cols, tokens.length);
                float total = Float.parseFloat(tokens[1]);
                float subs = Float.parseFloat(tokens[2]) + Float.parseFloat(tokens[3]);
                assertEquals(total, subs, 1e-2);
            }
        }
    }


    @Test
    public void testIncludeDuplicatesFlag() throws IOException {
        String bamURL = "http://www.broadinstitute.org/igvdata/BodyMap/hg18/Merged/HBM.adipose.bam.sorted.bam";
        String options = "d,q@chr1:153425249-153425249";
        int windowSize = 1;
        File wigFile = null;
        Genome genome = null;

        TestDataConsumer dc = new TestDataConsumer();

        CoverageCounter cc = new CoverageCounter(bamURL, dc, windowSize, 0, wigFile, genome, options);

        cc.parse();

        String totalCount = dc.attributes.get("totalCount");

        assertEquals("22", totalCount);

    }


    static class TestDataConsumer implements DataConsumer {

        Map<String, String> attributes = new HashMap<String, String>();

        public void setType(String type) {

        }

        public void addData(String chr, int start, int end, float[] data, String name) {

        }

        public void parsingComplete() {

        }

        public void setTrackParameters(TrackType trackType, String trackLine, String[] trackNames) {

        }

        public void setSortTolerance(int tolerance) {

        }

        public void setAttribute(String key, String value) {
            attributes.put(key, value);
        }
    }
}
