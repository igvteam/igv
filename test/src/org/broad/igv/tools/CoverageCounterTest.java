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

package org.broad.igv.tools;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.tools.parsers.DataConsumer;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.TestUtils;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import static junit.framework.Assert.assertEquals;


public class CoverageCounterTest extends AbstractHeadlessTest {

    static PreferenceManager preferenceManager;

    @BeforeClass
    public static void setUpClass() throws Exception {
        AbstractHeadlessTest.setUpClass();
        preferenceManager = PreferenceManager.getInstance();
    }

    /**
     * Test the "mapping quality" flag.  Also indirectly tests the query parameters.
     */
    @Test
    public void testMappingQualityFlag() throws IOException {
        String bamURL = TestUtils.LARGE_DATA_DIR + "HG00171.hg18.bam";
        String queryString = "chr1:152522155-152522155";
        int minMapQuality = 40;
        File wigFile = new File(TestUtils.DATA_DIR + "out/testMapQual.wig");
        int windowSize = 1;

        TestDataConsumer dc = new TestDataConsumer();

        CoverageCounter cc = new CoverageCounter(bamURL, dc, windowSize, 0, wigFile, genome, queryString, minMapQuality, 0);

        cc.parse();

        String totalCount = dc.attributes.get("totalCount");

        assertEquals("7", totalCount);
    }

    @Ignore    // The test file no longer exists
    @Test
    public void testPairFlag() throws Exception{
        String bamURL = "http://www.broadinstitute.org/igvdata/1KG/pilot2Bams/NA12878.SLX.bam";
        String queryString = "2:1000-1100";
        File wigFile = new File(TestUtils.DATA_DIR + "out/testPair.wig");
        int windowSize = 1;

        TestDataConsumer dc = new TestDataConsumer();

        CoverageCounter cc = new CoverageCounter(bamURL, dc, windowSize, 0, wigFile, genome, queryString, 0, CoverageCounter.PAIRED_COVERAGE);

        cc.parse();

        //Have manually checked these regions and verified that there is 1 pair
        //in the 1000-1100 region, and 4 at location 851
        //with
        for(TestData td: dc.testDatas){
            assertEquals(1.0f, td.data[0]);
        }


        queryString = "2:851";
        dc = new TestDataConsumer();

        cc = new CoverageCounter(bamURL, dc, windowSize, 0, wigFile, genome, queryString, 0, CoverageCounter.PAIRED_COVERAGE);

        cc.parse();

        for(TestData td: dc.testDatas){
            assertEquals(4.0f, td.data[0]);
        }

    }

    /*
    Test whether we count the number of features correctly.
     */
    @Test
    public void testCountStrand() throws Exception {
        String ifile = TestUtils.DATA_DIR + "bed/Unigene.sample.sorted.bed";


        int expTot = 71;
        int windowSize = 25;

        File wigFile = null;

        int[] countFlags = new int[]{0, CoverageCounter.STRANDS_BY_READ, CoverageCounter.STRANDS_BY_FIRST_IN_PAIR,
                CoverageCounter.BASES,
                CoverageCounter.INCLUDE_DUPS, CoverageCounter.BASES + CoverageCounter.STRANDS_BY_READ,
                CoverageCounter.BASES + CoverageCounter.STRANDS_BY_FIRST_IN_PAIR};
        //We specifically do not test this, because it's unreliable. Some alignments default to assuming they
        //are the first in pair, but let secondinpair be none
        //countFlags = new int[]{CoverageCounter.STRANDS_BY_SECOND_IN_PAIR};
        int[] expectedTotal = new int[]{expTot, expTot, expTot, expTot, expTot, expTot, expTot, expTot};

        for (int ii = 0; ii < countFlags.length; ii++) {
            TestDataConsumer dc = new TestDataConsumer();
            CoverageCounter cc = new CoverageCounter(ifile, dc, windowSize, 0, wigFile, genome, null, 0, countFlags[ii]);
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
        String ifile = TestUtils.DATA_DIR + "bed/Unigene.sample.sorted.bed";
        int[] windowSizes = new int[]{10, 50, 101, 500, 999};

        File wigFile = null;//new File(TestUtils.DATA_DIR + "out", "testStrandsConsistent.wig");
        //Test that when we run the process twice, with separate and totalled strands, the results add
        //up properly
        int[] strandOptions = new int[]{0, CoverageCounter.STRANDS_BY_READ, CoverageCounter.STRANDS_BY_FIRST_IN_PAIR, CoverageCounter.BASES};
        int[] expected_cols = new int[]{1, 2, 2, 5};
        TestDataConsumer[] tdcs = new TestDataConsumer[expected_cols.length];

        for (int ii = 0; ii < windowSizes.length; ii++) {
            for (int so = 0; so < strandOptions.length; so++) {
                TestDataConsumer dc = new TestDataConsumer();
                CoverageCounter cc = new CoverageCounter(ifile, dc, windowSizes[ii], 0, wigFile, genome, null, 0, strandOptions[so]);
                cc.parse();

                for (TestData tdata : dc.testDatas) {
                    float[] numbers = tdata.data;
                    assertEquals(expected_cols[so], numbers.length);
                }
                tdcs[so] = dc;
            }

            TestDataConsumer total_dc = tdcs[0];
            assertEquals(total_dc.testDatas.size(), tdcs[1].testDatas.size());
            for (int opts = 1; opts < strandOptions.length; opts++) {
                TestDataConsumer tdc = tdcs[opts];
                for (int row = 0; row < total_dc.testDatas.size(); row++) {
                    TestData td = tdc.testDatas.get(row);
                    float act_sum = 0;
                    for (float f : td.data) {
                        act_sum += f;
                    }
                    assertEquals(total_dc.testDatas.get(row).data[0], act_sum, 1e-2);
                }
            }
        }
    }

    @Test
    public void testCountBases() throws Exception {
        String ifile = TestUtils.DATA_DIR + "sam/NA12878.muc1.test.sam";
        int expected_cols = 10;

        File wigFile = new File(TestUtils.DATA_DIR + "out", "testCountBases.wig");
        int windowSize = 1;

        TestDataConsumer dc = new TestDataConsumer();
        int strandOptions = CoverageCounter.STRANDS_BY_READ + CoverageCounter.BASES;
        CoverageCounter cc = new CoverageCounter(ifile, dc, windowSize, 0, wigFile, genome, null, 0, strandOptions);
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

            if (tdata.start == check_startpos) {
                for (int ii = 0; ii < posvals.length; ii++) {
                    assertEquals(posvals[ii], numbers[ii], 1e-2);
                    assertEquals(negvals[ii], numbers[ii + keys.length], 1e-2);
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

    public void testColumnCounts() throws Exception {
        String ifile = TestUtils.DATA_DIR + "sam/NA12878.muc1.test.sam";

        File wigFile = null;
        Genome genome = this.genome;
        int[] windowSizes = new int[]{1, 50, 100, 500};
        //All possible combinations of STRAND_XXX flags
        int[] strandops = new int[2];
        strandops[0] = 0;
        strandops[1] = CoverageCounter.STRANDS_BY_READ;

//        String[] otherflags = new String[]{CoverageCounter.FIRST_IN_PAIR, CoverageCounter.BASES,
//                CoverageCounter.FIRST_IN_PAIR + CoverageCounter.BASES};

//        for (String so : strandops) {
//            for (String of : otherflags) {
//                int expectedcols = so.contains(CoverageCounter.STRAND_SEPARATE) ? 2 : 1;
//                if (of.contains(CoverageCounter.BASES)) {
//                    expectedcols *= 5;
//                }
//
//                String strandOptions = so + of;
//                for (int windowSize : windowSizes) {
//                    TestDataConsumer dc = new TestDataConsumer();
//                    CoverageCounter cc = new CoverageCounter(ifile, dc, windowSize, 0, wigFile, genome, "sc=" + strandOptions);
//                    cc.parse();
//
//                    assertEquals(expectedcols, dc.testDatas.get(0).data.length);
//                }
//            }
//        }

    }

    @Test
    public void testIncludeDuplicatesFlag() throws IOException {
        String bamURL = "http://www.broadinstitute.org/igvdata/BodyMap/hg18/Merged/HBM.adipose.bam.sorted.bam";
        int options = CoverageCounter.INCLUDE_DUPS;
        String queryString = "chr1:153425249-153425249";
        int windowSize = 1;
        File wigFile = null;
        Genome genome = null;

        TestDataConsumer dc = new TestDataConsumer();

        CoverageCounter cc = new CoverageCounter(bamURL, dc, windowSize, 0, wigFile, genome, queryString, 0, options);

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

        public void setTrackParameters(TrackType trackType, String trackLine, String[] trackNames, boolean b) {
            //To change body of implemented methods use File | Settings | File Templates.
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
