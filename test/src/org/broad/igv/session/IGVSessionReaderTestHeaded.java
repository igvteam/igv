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

package org.broad.igv.session;

import org.broad.igv.cli_plugin.PluginSpecReader;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.sam.CoverageTrack;
import org.broad.igv.track.DataSourceTrack;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.AbstractHeadedTest;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.TestUtils;
import org.broad.igv.variant.VariantTrack;
import org.broad.tribble.Feature;
import org.junit.Assume;
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
        String path = TestUtils.DATA_DIR + "sessions/coverage_snpThreshold.xml";
        rewriteRestoreSession(path);

        //We should have 1 coverage track
        CoverageTrack track = (CoverageTrack) IGV.getInstance().getAllTracks().get(0);

        assertTrue(track.isShowReference());
        assertEquals(0.5, track.getSnpThreshold(), 1e-5);
        assertTrue(track.isAutoScale());
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
        assertEquals(CoverageTrack.DEFAULT_AUTOSCALE, covTrack.isAutoScale());
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



    /**
     * Test loading a session with an analysis track
     * @throws Exception
     */
    @Test
    public void testLoadPluginSession() throws Exception{
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
        FeatureTrack matchTrack = (FeatureTrack) tracks.get(2);

        String queryChr = "chr8";
        int queryStart = 40823007;
        int queryEnd = 40863995;

        List<Feature> features = matchTrack.getFeatures(queryChr, queryStart, queryEnd);
        assertEquals(2, features.size());

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
