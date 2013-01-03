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
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

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
     * IGV-1421 - error reading a version 1 file (error thrown in header)
     */
    @Test
    public void test_IGV1417() {
        String tdfFile = "http://www.broadinstitute.org/igvdata/annotations/hg18/conservation/pi.ewig.tdf";
        TDFReader reader = TDFReader.getReader(tdfFile);
        assertEquals(1, reader.getVersion());
    }

    /**
     * IGV-1421 - error reading a version 2 file (error thrown in header)
     */
    @Test
    public void test_IGV1421() {
        String tdfFile = "http://www.broadinstitute.org/igvdata/1KG/freeze5/NA12878.pilot2.454.bam.tdf";
        TDFReader reader = TDFReader.getReader(tdfFile);
        assertEquals(2, reader.getVersion());
    }



    String[] dm3posChromos = new String[]{"chr2RHet", "chr4", "chrU"};
    String[] dm3emptyChromos = new String[]{"chrUextra"};

    @Test
    public void testChrAlldm3_v3() throws Exception{
        //TODO Put in test dir
        String genPath = "http://igvdata.broadinstitute.org/genomes/dm3.genome"; //"dm3";

        String wigPath = TestUtils.DATA_DIR + "wig/dm3_var_sample.wig";

        //TDF file generated from wiggle, using IGV 2.1.30 (tdf version 3)
        String tdf3Path = TestUtils.DATA_DIR + "tdf/dm3_var_sample.wig.v2.1.30.tdf";

        tstCHR_ALL(genPath, wigPath, tdf3Path, false, dm3posChromos, dm3emptyChromos);
    }

    @Test
    public void testChrAlldm3_v4() throws Exception{
        String genPath = "http://igvdata.broadinstitute.org/genomes/dm3.genome";

        String wigPath = TestUtils.DATA_DIR + "wig/dm3_var_sample.wig";
        String tdf4Path = TestUtils.DATA_DIR + "tdf/dm3_var_sample.wig.v2.2.1.tdf";

        tstCHR_ALL(genPath, wigPath, tdf4Path, true, dm3posChromos, dm3emptyChromos);
    }

    String[] hg18posChromos =  new String[]{"chr6", "chr7", "chr10", "chr11", "chr12"};
    String[] hg18emptyChromos = new String[]{"chr1", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrY"};
    @Test
    public void testChrAllhg18_v3() throws Exception{
        String genPath = TestUtils.DATA_DIR + "genomes/hg18.unittest.genome";

        String wigPath = TestUtils.DATA_DIR + "wig/hg18_var_sample.wig";

        //TDF file generated from wiggle, using IGV 2.1.30 (tdf version 3)
        String tdf3Path = TestUtils.DATA_DIR + "tdf/hg18_var_sample.wig.v2.1.30.tdf";

        //Chromosome order for this genome didn't change, at least for the relevant chromosomes
        tstCHR_ALL(genPath, wigPath, tdf3Path, true, hg18posChromos, hg18emptyChromos);
    }

    @Test
    public void testChrAllhg18_v4() throws Exception{
        String genPath = TestUtils.DATA_DIR + "genomes/hg18.unittest.genome";

        String wigPath = TestUtils.DATA_DIR + "wig/hg18_var_sample.wig";
        String tdf4Path = TestUtils.DATA_DIR + "tdf/hg18_var_sample.wig.v2.2.1.tdf";

        tstCHR_ALL(genPath, wigPath, tdf4Path, true, hg18posChromos, hg18emptyChromos);
    }

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
        Genome genome = GenomeManager.getInstance().loadGenome(genPath, null);

        WiggleDataset ds = (new WiggleParser(new ResourceLocator(wigPath), genome)).parse();
        DatasetDataSource wigSource = new DatasetDataSource(wigPath, ds, genome);


        TDFReader tdfReader = TDFReader.getReader(tdfPath);
        TDFDataSource tdfSource = new TDFDataSource(tdfReader, 0, tdfPath, genome);

        //We can't test for exact equality, so we just look for errors
        List<LocusScore> wigScores = wigSource.getSummaryScoresForRange(Globals.CHR_ALL, -1, -1, 0);
        List<LocusScore> tdfScores = tdfSource.getSummaryScoresForRange(Globals.CHR_ALL, -1, -1, 0);

        assertEquals(expHaveChrAll, tdfSource.isChrAllValid());

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

        assertTrue(posChecked > 0);
        //System.out.println("# Checked for positive values: " + posChecked);

        for(String chromo: emptyChromos){
            long range = genome.getChromosome(chromo).getLength() / 1000;
            long fudge = Math.max(500, range / 100);
            range -= 2 * fudge;

            long cMin = (genome.getCumulativeOffset(chromo) / 1000) + fudge;
            long cMax = cMin + range;

            for(LocusScore wLocus: wigScores){
                //Just checking our assumption
                assert wLocus.getStart() < cMin || wLocus.getEnd() > cMax;
            }

            for(LocusScore tdfScore: tdfScores){
                assertFalse("Found data where none should exist at " + chromo + ":" + tdfScore.getStart() + "-" + tdfScore.getEnd(), overlaps(cMin, cMax, tdfScore));
            }

        }
    }

}
