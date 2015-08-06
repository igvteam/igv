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

package org.broad.igv.session;

import com.google.common.base.Predicate;
import com.google.common.collect.Iterables;
import org.broad.igv.cli_plugin.AbstractPluginTest;
import org.broad.igv.cli_plugin.PluginSpecReader;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.gwas.GWASTrack;
import org.broad.igv.renderer.BarChartRenderer;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.sam.CoverageTrack;
import org.broad.igv.track.*;
import org.broad.igv.ui.AbstractHeadedTest;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.TestUtils;
import org.broad.igv.variant.VariantTrack;
import htsjdk.tribble.Feature;
import org.junit.Assume;
import org.junit.Ignore;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.*;

/** Test reading sessions, make sure attributes get parsed properly
 *
 * TODO USE LOCAL DATA FILES
 * User: jacob
 * Date: 2013-Jan-07
 */
public class IGVSessionReaderTestHeaded extends AbstractHeadedTest{

    @Test
    public void testLoadCoverageTrackSession() throws Exception{
        String path = TestUtils.DATA_DIR + "sessions/coverage_snpThreshold_NA06984.xml";
        rewriteRestoreSession(path);

        //We should have 1 coverage track
        CoverageTrack track = (CoverageTrack) IGV.getInstance().getAllTracks().get(0);

        assertTrue(track.isShowReference());
        assertEquals(0.5, track.getSnpThreshold(), 1e-5);
        assertTrue(track.getAutoScale());
    }

    @Test
    public void testLoadDataSourceTrack() throws Exception{
        String path = TestUtils.DATA_DIR + "sessions/tdf_session.xml";
        rewriteRestoreSession(path);

        DataSourceTrack track = (DataSourceTrack) IGV.getInstance().getAllTracks().get(0);
        assertEquals(true, TestUtils.runMethod(track, "getNormalize"));
        assertEquals("NA12878.SLX.egfr.sam.tdf", track.getName());
    }

    /**
     * Loads path, returns first alignment track
     * @param sessionPath
     * @return
     */
    private AlignmentTrack getAlignmentTrack(String sessionPath) throws Exception{
        rewriteRestoreSession(sessionPath);

        AlignmentTrack alTrack = null;
        for(Track tk: IGV.getInstance().getAllTracks()){
            if(tk instanceof AlignmentTrack){
                alTrack = (AlignmentTrack) tk;
            }
        }

        return alTrack;
    }

    @Test
    public void testLoadAlignmentTrackSession() throws Exception{
        String path = TestUtils.DATA_DIR + "sessions/HG00171_wlocus_v5.xml";
        AlignmentTrack alTrack = getAlignmentTrack(path);

        CoverageTrack covTrack = alTrack.getCoverageTrack();
        assertTrue(alTrack.isShowSpliceJunctions());

        assertEquals(CoverageTrack.DEFAULT_SHOW_REFERENCE, covTrack.isShowReference());
        assertEquals(CoverageTrack.DEFAULT_AUTOSCALE, covTrack.getAutoScale());
        assertEquals(0.1337f, covTrack.getSnpThreshold(), 1e-5);
    }

    @Test
    public void testLoadVCF_v4() throws Exception{
        String path = TestUtils.DATA_DIR + "sessions/vcf_session_v4.xml";
        tstLoadVCF(path);
    }

    @Test
    public void testLoadVCF_v5() throws Exception{
        String path = TestUtils.DATA_DIR + "sessions/vcf_session_v5.xml";
        tstLoadVCF(path);
    }

    private void tstLoadVCF(String sessionPath) throws Exception{
        rewriteRestoreSession(sessionPath);

        VariantTrack vTrack = (VariantTrack) IGV.getInstance().getAllTracks().get(0);

        assertEquals(3, vTrack.getAllSamples().size());
        assertEquals(8, vTrack.getSquishedHeight());
        assertEquals(VariantTrack.ColorMode.ALLELE, vTrack.getColorMode());
    }


    @Test
    public void testLoadRenderOptions_v4() throws Exception{
        String path = TestUtils.DATA_DIR + "sessions/HG00171_wlocus_v4.xml";
        tstLoadRenderOptions(path);
    }

    @Test
    public void testLoadRenderOptions_v5() throws Exception{
        String path = TestUtils.DATA_DIR + "sessions/HG00171_wlocus_v5.xml";
        tstLoadRenderOptions(path);
    }

    private void tstLoadRenderOptions(String path) throws Exception{
        AlignmentTrack alTrack = getAlignmentTrack(path);

        AlignmentTrack.RenderOptions ro = getRenderOptions(alTrack);

        assertEquals(AlignmentTrack.GroupOption.FIRST_OF_PAIR_STRAND, ro.getGroupByOption());
        assertEquals(AlignmentTrack.ColorOption.INSERT_SIZE, ro.getColorOption());
        assertEquals(55, ro.getMinInsertSize());
        assertEquals(1001, ro.getMaxInsertSize());
        assertEquals("", ro.getColorByTag());
    }

    private static class ContIdPredicate implements Predicate<Track>{

        private String[] contIds;

        ContIdPredicate(String contId){
            this.contIds = new String[]{contId};
        }

        ContIdPredicate(String[] contIds){
            this.contIds = contIds;
        }

        @Override
        public boolean apply(Track input) {
            boolean contains = true;
            for(String contId: contIds){
                contains &= input.getId().contains(contId);
                if(!contains) break;
            }
            return contains;
        }

    }

    @Test
    public void testLoadCombinedDataSourceSession() throws Exception{
        String sessionpath = TestUtils.DATA_DIR + "sessions/subtypes_wdiff.xml";
        rewriteRestoreSession(sessionpath);

        String combPart0 = "TRIBE_p_TCGAaffx_B1_2_GBM_Nsp_GenomeWideSNP_6_E11_155884";
        String combPart1 = "TRIGS_p_TCGAaffxB5_sty_GenomeWideSNP_6_D05_223156";

        DataTrack track0 = Iterables.find(IGV.getInstance().getDataTracks(), new ContIdPredicate(combPart0));
        DataTrack track1 = Iterables.find(IGV.getInstance().getDataTracks(), new ContIdPredicate(combPart1));

        DataSourceTrack combTrack = (DataSourceTrack) Iterables.find(IGV.getInstance().getDataTracks(), new ContIdPredicate(new String[]{combPart0, combPart1}));

        assertTrue(combTrack.getRenderer() instanceof BarChartRenderer);
        assertEquals("Difference", combTrack.getName());


        String chr = "chr1";
        int start = 0;
        int end = 1000;

        LocusScore sumScore0 = track0.getSummaryScores(chr, start, end, 0).get(0);
        LocusScore sumScore1 = track1.getSummaryScores(chr, start, end, 0).get(0);
        LocusScore sumScoreComb = combTrack.getSummaryScores(chr, start, end, 0).get(0);

        assertEquals(sumScore0.getStart(), sumScore1.getStart());
        assertEquals(sumScore0.getStart(), sumScoreComb.getStart());
        assertEquals(sumScore1.getEnd(), sumScoreComb.getEnd());

        assertEquals(sumScore0.getScore() - sumScore1.getScore(), sumScoreComb.getScore(), 1e-10);

    }


    /**
     * Test creating an analysis track, saving the session, and loading it
     */
    @Ignore("Not ready yet")
    @Test
    public void testSaveLoadPluginSession() throws Exception{
        String toolPath = "/usr/local/bin/bedtools";
        boolean haveTool = PluginSpecReader.isToolPathValid(toolPath);
        Assume.assumeTrue(haveTool);

        String sessionPath = TestUtils.DATA_DIR + "sessions/GSM_beds.xml";
        rewriteRestoreSession(sessionPath);

        String trackAname = "GSM1004654_100k.bed";
        String trackBname = "GSM1004654_10k.bed";

        //Generate bedtools subtraction track
        PluginSpecReader reader = PluginSpecReader.create(toolPath);
        PluginSpecReader.Tool tool = reader.getTools().get(0);
        PluginSpecReader.Command command = AbstractPluginTest.findCommandElementByName(tool, "Subtract");
    }

    @Test
    public void testLoadOverlaySession() throws Exception{
        String sessionpath = TestUtils.DATA_DIR + "sessions/hind_gistic_overlay.xml";
        rewriteRestoreSession(sessionpath);

        //Load the session, check things went well
        List<DataTrack> dataTracks = IGV.getInstance().getDataTracks();
        assertEquals(3, dataTracks.size());
        assertTrue(dataTracks.get(0) instanceof MergedTracks);
    }

    /**
     * Very basic test, just makes sure it loads without throwing an exception
     * @throws Exception
     */
    @Test
    public void testLoadBedtools_Intersect() throws Exception{
        String sessionPath = TestUtils.DATA_DIR + "sessions/GSM_bedtools_intersect.xml";
        String toolPath = "/usr/local/bin/bedtools";
        boolean haveTool = PluginSpecReader.isToolPathValid(toolPath);
        Assume.assumeTrue(haveTool);

        rewriteRestoreSession(sessionPath);
    }
    /**
     * Test loading a session with a bedtools analysis track
     * @throws Exception
     */
    @Test
    public void testLoadBedtools_Subtract() throws Exception{
        String toolPath = "/usr/local/bin/bedtools";
        boolean haveTool = PluginSpecReader.isToolPathValid(toolPath);
        Assume.assumeTrue(haveTool);

        String sessionPath = TestUtils.DATA_DIR + "sessions/GSM_bedtools_subtract.xml";

        rewriteRestoreSession(sessionPath);

        String trackAname = "GSM1004654_100k.bed";
        String trackBname = "GSM1004654_10k.bed";
        String analId = "BEDTools Remove/Subtract";

        FeatureTrack analTrack = null, trackA = null, trackB = null;
        for(Track track: IGV.getInstance().getAllTracks()){
            if(track.getId().equals(analId)){
                analTrack = (FeatureTrack) track;
            }else if(track.getName().equals(trackAname)){
                trackA = (FeatureTrack) track;
            }else if(track.getName().equals(trackBname)){
                trackB = (FeatureTrack) track;
            }
        }

        String chr = "chr2";
        int start = 177932002;
        int end = 180561093;

        List<Feature> aFeatures = trackA.getFeatures(chr, start, end);
        List<Feature> bFeatures = trackB.getFeatures(chr, start, end);
        List<Feature> analFeatures = analTrack.getFeatures(chr, start, end);

        //Fairly coarse check, these actually might not exactly be true
        //due to splitting of exons
        int checked = 0;
        for (Feature afeat : analFeatures) {
            if(afeat.getStart() < start || afeat.getStart() > end) continue;
            //This particular feature splits funny, it's not a bug, at least not with IGV
            if(afeat.getStart() == 178625608) continue;

            assertTrue(listContainsFeature(aFeatures, afeat));
            assertFalse(listContainsFeature(bFeatures, afeat));
            checked++;
        }
        assert checked > 0;
    }


    /**
     * Test loading a session which contains a FeatureSource for
     * sequence matching.
     * @throws Exception
     */
    @Test
    public void testLoadFeatureMatchSession() throws Exception{
        String sessionPath = TestUtils.DATA_DIR + "sessions/match_search_session.xml";
        IGV.getInstance().doRestoreSession(sessionPath, null, false);

        List<Track> tracks = IGV.getInstance().getAllTracks();

        //Positive strand
        FeatureTrack matchTrack = (FeatureTrack) tracks.get(2);

        String queryChr = "chr8";
        int queryStart = 40823007;
        int queryEnd = 40863995;

        List<Feature> features = matchTrack.getFeatures(queryChr, queryStart, queryEnd);
        assertEquals(2, features.size());

        //Negative strand
        matchTrack = (FeatureTrack) tracks.get(3);
        features = matchTrack.getFeatures(queryChr, queryStart, queryEnd);
        assertEquals(1, features.size());

    }

    @Test
    public void testLoadNoGeneTrack() throws Exception{

        assertTrue(IGV.getInstance().hasGeneTrack());

        String sessionPath = TestUtils.DATA_DIR + "sessions/no_gene_track.xml";
        IGV.getInstance().doRestoreSession(sessionPath, null, false);

        assertFalse("Gene track exists but it shouldn't", IGV.getInstance().hasGeneTrack());

        sessionPath = TestUtils.DATA_DIR + "sessions/has_gene_track.xml";
        IGV.getInstance().doRestoreSession(sessionPath, null, false);

        assertTrue("No gene track found but one should exist", IGV.getInstance().hasGeneTrack());
    }


    String EcoliFastaSessionPath = TestUtils.DATA_DIR + "sessions/ecoli_out_seqonly.xml";

    @Test
    public void testLoadSequenceTrackOnlyDiffGenome() throws Exception{
        rewriteRestoreSession(EcoliFastaSessionPath);

        assertFalse("Gene track exists but it shouldn't", IGV.getInstance().hasGeneTrack());

        assertTrue("Sequence track doesn't exist but it should", IGV.getInstance().hasSequenceTrack());
        assertNotNull("Sequence track doesn't exist but it should", IGV.getInstance().getSequenceTrack());
    }

    @Test
    public void testLoadSequenceTrackOnlySameGenome() throws Exception{
        String genomePath = TestUtils.DATA_DIR + "fasta/ecoli_out.padded.fasta";
        GenomeManager.getInstance().loadGenome(genomePath, null);
        genome = GenomeManager.getInstance().getCurrentGenome();
        IGV.getInstance().removeTracks(IGV.getInstance().getAllTracks());

        assertFalse("Sequence Track exists but it shouldn't", IGV.getInstance().hasSequenceTrack());

        String sessionPath = TestUtils.DATA_DIR + "sessions/ecoli_out_seqonly.xml";
        rewriteRestoreSession(sessionPath);

        assertFalse("Gene track exists but it shouldn't", IGV.getInstance().hasGeneTrack());

        assertTrue("Sequence track doesn't exist but it should", IGV.getInstance().hasSequenceTrack());
        assertNotNull("Sequence track doesn't exist but it should", IGV.getInstance().getSequenceTrack());
    }


    /**
     * Test loading a session with spaces in the path and special
     * characters in the file names
     * @throws Exception
     */
    @Test
    public void testLoadSessionSpecialChars() throws Exception {
        int expTrackCount = 5;
        String path = TestUtils.DATA_DIR + "sessions/special_chars.xml";
        rewriteRestoreSession(path);

        assertEquals(expTrackCount, IGV.getInstance().getVisibleTrackCount());

    }

    @Test
    public void testLoadSessionWithoutLoadedGenome() throws Exception{
        GenomeManager.getInstance().setCurrentGenome(null);
        assertNull(GenomeManager.getInstance().getGenomeId());

        rewriteRestoreSession(EcoliFastaSessionPath);

        assertNotNull(GenomeManager.getInstance().getGenomeId());


    }

    @Test
    public void testGWAS() throws Exception{
        String sessionPath = TestUtils.DATA_DIR + "sessions/smallp_session.xml";
        rewriteRestoreSession(sessionPath);
        boolean haveGWAS = false;
        for(Track track: IGV.getInstance().getAllTracks()){
            if(track instanceof GWASTrack){
                haveGWAS = true;
                break;
            }
        }
        assertTrue("No GWAS track loaded", haveGWAS);
        assertEquals(3, IGV.getInstance().getVisibleTrackCount());
    }

    private boolean listContainsFeature(List<Feature> featureList, Feature feature){
        for (Feature listFeature : featureList) {
            if ((feature.getChr().equals(listFeature.getChr())) &&
                    (feature.getStart() == listFeature.getStart()) &&
                    (feature.getEnd() == listFeature.getEnd())) {
                return true;
            }
        }
        return false;
    }

    /**
     * Get render options of the track, reflectively so we don't need to alter access controls
     * @param track
     * @return
     * @throws Exception
     */
    private AlignmentTrack.RenderOptions getRenderOptions(Track track) throws Exception{
        Object result = TestUtils.getField(track, "renderOptions");
        return (AlignmentTrack.RenderOptions) result;
    }
}
