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

package org.broad.igv.tdf;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.Globals;
import org.broad.igv.data.DatasetDataSource;
import org.broad.igv.data.WiggleDataset;
import org.broad.igv.data.WiggleParser;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Assume;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import java.io.IOException;
import java.util.List;

import static org.junit.Assert.*;


/**
 * @author jrobinso
 * @date Jul 28, 2010
 */
public class TDFRegressionTests extends AbstractHeadlessTest{



    @Rule
    public TestRule testTimeout = new Timeout((int) 30e6);

     /**
     * IGV-1417 and/or IGV-1421 - error reading a version 1 file (error thrown in header)
     */
    @Test
    public void test_IGV1417() {
        String tdfFile = "http://www.broadinstitute.org/igvdata/annotations/hg18/conservation/pi.ewig.tdf";
        TDFReader reader = TDFReader.getReader(tdfFile);
        assertEquals(1, reader.getVersion());
    }

    @Test
    public void test_v3() {
        String tdfFile = "http://www.broadinstitute.org/igvdata/test/tdf/NA12878.pilot2.454.bam.tdf";
        TDFReader reader = TDFReader.getReader(tdfFile);
        assertEquals(3, reader.getVersion());
    }

    String[] dm3posChromos = new String[]{"chr2RHet", "chr4", "chrU"};
    String[] dm3emptyChromos = new String[]{"chrUextra"};

//    NOTE:  V3 files for genomes other than hg18, hg19, mm8, or mm9 always fail by definition, so there is nothing to test
//    @Test
//    public void testChrAlldm3_v3() throws Exception{
//        //TODO Put in test dir
//        String genPath = "http://igvdata.broadinstitute.org/genomes/dm3.genome"; //"dm3";
//
//        String wigPath = TestUtils.DATA_DIR + "wig/dm3_var_sample.wig";
//
//        //TDF file generated from wiggle, using IGV 2.1.30 (tdf version 3)
//        String tdf3Path = TestUtils.DATA_DIR + "tdf/dm3_var_sample.wig.v2.1.30.tdf";
//
//        tstCHR_ALL(genPath, wigPath, tdf3Path, false, dm3posChromos, dm3emptyChromos);
//    }

    @Test
    public void testChrAlldm3_v4() throws Exception{
        String genPath = "http://igvdata.broadinstitute.org/genomes/dm3.genome";

        String wigPath = TestUtils.DATA_DIR + "wig/dm3_var_sample.wig";
        String tdf4Path = TestUtils.DATA_DIR + "tdf/dm3_var_sample.wig.v2.2.1.tdf";

        // The chr order test should fail as order for dm3  has been changed with release v2.3
        tstCHR_ALL(genPath, wigPath, tdf4Path, false, dm3posChromos, dm3emptyChromos);
    }

    String[] hg18posChromos =  new String[]{"chr6", "chr7", "chr10", "chr11", "chr12"};
    String[] hg18emptyChromos = new String[]{"chr1", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrY"};

    @Test
    public void testChrAllhg18_v3() throws Exception{
        String genPath = TestUtils.DATA_DIR + "genomes/hg18.unittest.genome";

        String wigPath = TestUtils.DATA_DIR + "wig/hg18_var_sample.wig";

        //TDF file generated from wiggle, using IGV 2.1.30 (tdf version 3)
        String tdf3Path = TestUtils.DATA_DIR + "tdf/hg18_var_sample.wig.v2.1.30.tdf";

        tstCHR_ALL(genPath, wigPath, tdf3Path, false, hg18posChromos, hg18emptyChromos);
    }
//
//    @Test
//    public void testChrAllhg18_v4() throws Exception{
//        String genPath = TestUtils.DATA_DIR + "genomes/hg18.unittest.genome";
//
//        String wigPath = TestUtils.DATA_DIR + "wig/hg18_var_sample.wig";
//        String tdf4Path = TestUtils.DATA_DIR + "tdf/hg18_var_sample.wig.v2.2.1.tdf";
//
//        tstCHR_ALL(genPath, wigPath, tdf4Path, true, hg18posChromos, hg18emptyChromos);
//    }

    private boolean overlaps(long start, long end, LocusScore score){
        return score.getStart() >= start && score.getStart() < end ||
                score.getStart() <= end && score.getEnd() > start;
    }

    /**
     * Test that our whole genome view returns valid data, or none at all
     *
     * @param genPath Path to genome used for tdf
     * @param wigPath wig file used as input to tdf TODO Any toTDF file
     * @param tdfPath path to outputted tdf file
     * @param expHaveChrAll Whether the CHR_ALL call to the datasource is expected to be valid
     * @param posChromos Chromosomes for which all values are positive
     * @param emptyChromos Chromosomes which have no data
     * @throws Exception
     */
    public void tstCHR_ALL(String genPath, String wigPath, String tdfPath, boolean expHaveChrAll, String[] posChromos, String[] emptyChromos) throws Exception{
        Genome genome = null;
        try {
            genome = GenomeManager.getInstance().loadGenome(genPath, null);
        } catch (IOException e) {
            e.printStackTrace();
        }
        Assume.assumeTrue(genome != null);


        WiggleDataset ds = (new WiggleParser(new ResourceLocator(wigPath), genome)).parse();
        DatasetDataSource wigSource = new DatasetDataSource(wigPath, ds, genome);


        TDFReader tdfReader = TDFReader.getReader(tdfPath);
        TDFDataSource tdfSource = new TDFDataSource(tdfReader, 0, tdfPath, genome);

        //We can't test for exact equality, so we just look for errors
        List<LocusScore> wigScores = wigSource.getSummaryScoresForRange(Globals.CHR_ALL, -1, -1, 0);
        List<LocusScore> tdfScores = tdfSource.getSummaryScoresForRange(Globals.CHR_ALL, -1, -1, 0);

        assertEquals(expHaveChrAll, tdfSource.isChrOrderValid());

        if(!expHaveChrAll){
            //Ideally we would recalculate the data, but returning nothing
            //is preferable to returning incorrect data
            if(tdfScores.size() == 0){
                return;
            }
        }


        //TODO This test could be more efficient, and we could be more precise
        //about chromosome boundaries
        int posChecked = 0;
        for(String chromo: posChromos){
            //We allow a little slop around chromo boundaries
            long range = genome.getChromosome(chromo).getLength() / 1000;
            long fudge = Math.max(500, range / 100);
            range -= 2 * fudge;

            long cMin = genome.getCumulativeOffset(chromo) / 1000 + fudge;
            long cMax = cMin + range;

            for(int ii= 0; ii < wigScores.size(); ii++){
                LocusScore wLocus = wigScores.get(ii);
                if(wLocus.getStart() >= cMin && wLocus.getEnd() <= cMax){
                    float wScore = wigScores.get(ii).getScore();
                    //Just checking our assumption
                    assert wScore > 0;
                }

                LocusScore tLocus = tdfScores.get(ii);
                if(overlaps(cMin, cMax, tLocus)){
                    assertTrue("Found negative value in " + chromo + " at " + tLocus.getStart(), tLocus.getScore() >= 0);
                    posChecked++;
                }
            }

        }

        //assertTrue(posChecked > 0);
        System.out.println("# Checked for positive values: " + posChecked);

        for(String chromo: emptyChromos){
            long range = genome.getChromosome(chromo).getLength() / 1000;

            long fudge = Math.max(600, range / 100);
            range -= 2 * fudge;

            long cMin = (genome.getCumulativeOffset(chromo) / 1000) + fudge;
            long cMax = cMin + range;

            for(LocusScore wLocus: wigScores){
                //Just checking our assumption
                assert wLocus.getStart() < cMin || wLocus.getEnd() > cMax;
            }

            for(LocusScore tdfScore: tdfScores){
                boolean hasData = overlaps(cMin, cMax, tdfScore);
                hasData &= tdfScore.getScore() != 0.0f;
                assertFalse("Found data where none should exist at " + chromo + ":" + tdfScore.getStart() + "-" + tdfScore.getEnd(), hasData);
            }

        }
    }

}
