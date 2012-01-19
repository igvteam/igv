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
import org.broad.igv.util.Utilities;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
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

        File wigFile = null;
        Genome genome = this.genome;


        String[] options = new String[]{"sc=0x01", "sc=0x04", "sc=0x08", "sc=0xC", "sc=0x1C", "sc=0x05", "sc=0x02"};
        int[] expectedTotal = new int[]{expTot, expTot, expTot, expTot, expTot, expTot, expTot};

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
        int expected_cols = 3;

        File wigFile = null;//new File(TestUtils.DATA_DIR + "/out", "testStrandsConsistent.wig");
        Genome genome = this.genome;
        int strandOptions = CoverageCounter.STRAND_TOTAL + CoverageCounter.STRAND_ZERO + CoverageCounter.STRAND_ONE;


        for (int ii = 0; ii < windowSizes.length; ii++) {
            TestDataConsumer dc = new TestDataConsumer();
            CoverageCounter cc = new CoverageCounter(ifile, dc, windowSizes[ii], 0, wigFile, genome, "sc=" + strandOptions);
            cc.parse();

            for (TestData tdata : dc.testDatas) {
                float[] numbers = tdata.data;
                assertEquals(expected_cols, numbers.length);
                float total = numbers[0];
                float subs = numbers[1] + numbers[2];
                assertEquals(total, subs, 1e-2);
            }
        }
    }

    @Test
    public void testCountBases() throws Exception {
        String ifile = TestUtils.DATA_DIR + "/sam/NA12878.muc1.test.sam";
        int expected_cols = 18;

        File wigFile = new File(TestUtils.DATA_DIR + "/out", "testCountBases.wig");
        Genome genome = this.genome;
        int windowSize = 1;

        TestDataConsumer dc = new TestDataConsumer();
        int strandOptions = CoverageCounter.STRAND_TOTAL + CoverageCounter.STRAND_ZERO + CoverageCounter.STRAND_ONE + CoverageCounter.BASES;
        CoverageCounter cc = new CoverageCounter(ifile, dc, windowSize, 0, wigFile, genome, "sc=" + strandOptions);
        cc.parse();


        int check_startpos = 153426135 - 1;
        Map<Byte, Integer> posCounts = new HashMap<Byte, Integer>();
        Map<Byte, Integer> negCounts = new HashMap<Byte, Integer>();
        byte[] keys = new byte[]{'A', 'C', 'G', 'T', 'N'};
        int[] posvals = new int[]{9, 0, 0, 2, 0};
        int[] negvals = new int[]{16, 0, 0, 0, 0};
        for (int ii = 0; ii < keys.length; ii++) {
            posCounts.put(keys[ii], posvals[ii]);
            negCounts.put(keys[ii], negvals[ii]);
        }

        for (TestData tdata : dc.testDatas) {

            float[] numbers = tdata.data;
            assertEquals(expected_cols, numbers.length);
            float total = numbers[0];
            float subs = 0;
            for (int ii = 1; ii <= 2; ii++) {
                //Check that they add up total in each row
                //just looking at strands
                subs += numbers[ii];
            }
            assertEquals(total, subs, 1e-2);

            //Count up bases for each strand, test
            //that they add up right
            subs = 0;
            for (int ii = 3; ii < 3 + 5; ii++) {
                assertEquals(numbers[ii], numbers[ii + 5] + numbers[ii + 2 * 5], 1e-2);
            }

            if (tdata.start == check_startpos) {
                for (int ii = 0; ii < posvals.length; ii++) {
                    assertEquals(posvals[ii], numbers[ii + 3 + 5], 1e-2);
                    assertEquals(negvals[ii], numbers[ii + 3 + 2 * 5], 1e-2);
                }
            }
        }

    }

    /**
     * Test different strand options, just count output columns
     * and make sure we get the right number
     *
     * @throws Exception
     */
    @Test
    public void testColumnCounts() throws Exception {
        String ifile = TestUtils.DATA_DIR + "/sam/NA12878.muc1.test.sam";

        File wigFile = null;
        Genome genome = this.genome;
        int[] windowSizes = new int[]{1, 50, 100, 500};
        //All possible combinations of STRAND_XXX flags
        int[] strandops = new int[7];
        strandops[0] = CoverageCounter.STRAND_ZERO;
        strandops[1] = CoverageCounter.STRAND_ONE;
        strandops[2] = CoverageCounter.STRAND_TOTAL;
        strandops[3] = strandops[0] + strandops[1];
        strandops[4] = strandops[0] + strandops[2];
        strandops[5] = strandops[1] + strandops[2];
        strandops[6] = strandops[3] + strandops[2];

        int[] otherflags = new int[]{CoverageCounter.FIRST_IN_PAIR, CoverageCounter.BASES,
                CoverageCounter.FIRST_IN_PAIR + CoverageCounter.BASES};

        for (int so : strandops) {
            for (int of : otherflags) {
                int expectedcols = Utilities.countFlags(so);
                if ((of & CoverageCounter.BASES) > 0) {
                    expectedcols *= 6;
                }

                int strandOptions = so + of;
                for (int windowSize : windowSizes) {
                    TestDataConsumer dc = new TestDataConsumer();
                    CoverageCounter cc = new CoverageCounter(ifile, dc, windowSize, 0, wigFile, genome, "sc=" + strandOptions);
                    cc.parse();

                    assertEquals(expectedcols, dc.testDatas.get(0).data.length);
                }
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
        public ArrayList<TestData> testDatas = new ArrayList<TestData>();

        public void setType(String type) {

        }

        /**
         * Just write all the data a containment object
         *
         * @param chr
         * @param start
         * @param end
         * @param data
         * @param name
         */
        public void addData(String chr, int start, int end, float[] data, String name) {
            TestData newData = new TestData(chr, start, end, data.clone(), name);
            this.testDatas.add(newData);
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

    private static class TestData {
        public String chr;
        public int start;
        public int end;
        public float[] data;
        public String name;

        public TestData(String chr, int start, int end, float[] data, String name) {
            this.chr = chr;
            this.start = start;
            this.end = end;
            this.data = data;
            this.name = name;
        }
    }
}
