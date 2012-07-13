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

package org.broad.igv.batch;

import org.broad.igv.AbstractHeadedTest;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.IGVTestHeadless;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.util.TestUtils;
import org.junit.*;

import java.io.File;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import static junit.framework.Assert.*;

/**
 * User: jacob
 * Date: 2012/03/21
 */
public class CommandExecutorTest extends AbstractHeadedTest{

    CommandExecutor exec = new CommandExecutor();
    private final String snapshotDir = TestUtils.DATA_DIR + "out/";

    @Before
    public void setUp() throws Exception {
        super.setUp();
        igv.getSession().clearRegionsOfInterest();
        exec.setSnapshotDirectory(snapshotDir);
    }

    @After
    public void tearDown() throws Exception {
        igv.removeTracks(igv.getAllTracks());
    }

    @Test
    public void testRegion() throws Exception {

        String regionStr = "chr1:50-1000";
        exec.execute("region " + regionStr);

        Collection<RegionOfInterest> regions = IGV.getInstance().getSession().getAllRegionsOfInterest();
        assertEquals(1, regions.size());

        ArrayList<RegionOfInterest> regionsAL = (ArrayList<RegionOfInterest>) regions;
        RegionOfInterest region = regionsAL.get(0);
        assertEquals("chr1", region.getChr());
        assertEquals(49, region.getStart());
        assertEquals(1000, region.getEnd());

    }

    @Test
    public void testSetMaxDepth() throws Exception {
        setCheckMaxDepth(5);
        setCheckMaxDepth(50);
    }

    private void setCheckMaxDepth(int maxDepth) {
        String res = exec.execute("maxDepth " + maxDepth);
        assertFalse(res.contains("ERROR"));
        int newMaxDepth = PreferenceManager.getInstance().getAsInt(PreferenceManager.SAM_SAMPLING_COUNT);
        assertEquals(maxDepth, newMaxDepth);

    }

    private List<AlignmentTrack> getAlTracks(List<Track> tracks) {
        List<AlignmentTrack> out = new ArrayList<AlignmentTrack>(tracks.size());
        for (Track t : tracks) {
            try {
                out.add((AlignmentTrack) t);
            } catch (ClassCastException e) {
                continue;
            }
        }

        return out;
    }

    @Test
    public void testSortByRegionScoreType() throws Exception {
        String sessionPath = TestUtils.DATA_DIR + "sessions/BRCA_loh2.xml";
        TestUtils.loadSession(igv, sessionPath);
        Collection<RegionOfInterest> rois = igv.getSession().getAllRegionsOfInterest();

        //assertEquals(1, rois.size());

        List<Track> tracks;
        int count = 0;
        for (RegionOfInterest roi : rois) {
            for (RegionScoreType type : RegionScoreType.values()) {
                igv.sortAllTracksByAttributes(new String[]{"NAME"}, new boolean[]{false});
                String typeStr = type.toString().toUpperCase();
                if (count % 2 == 0) {
                    typeStr = typeStr.toLowerCase();
                }
                String resp = exec.execute("sort " + typeStr + " " + roi.getLocusString());
                assertEquals("OK", resp);

                tracks = igv.getAllTracks();
                IGVTestHeadless.checkIsSorted(tracks, roi, type, FrameManager.getDefaultFrame().getZoom());
                count++;
            }
        }


    }

    private final String outFileBase = "testSnap";

    @Test
    public void testSnapShotPng() throws Exception{
        String outFileName =  outFileBase + ".Png";
        tstSnapshot(outFileName);
    }

    @Test
    public void testSnapShotJpeg() throws Exception{
        tstSnapshot(outFileBase + ".jpeg");
    }

    @Test
    public void testSnapShotJpg() throws Exception{
        tstSnapshot(outFileBase + ".jpg");
    }

    @Test
    public void testSnapShotSvg() throws Exception{
        String outFileName =  outFileBase + ".svG";
        tstSnapshot(outFileName);
    }

    @Test
    public void testSnapShotFails() throws Exception{
        String[] exts = new String[]{"abc", "svt", "pnq"};
        for(String ext: exts){
            String outFileName =  outFileBase + "." + ext;
            tstSnapshot(outFileName, false);
        }
    }


    public void tstSnapshot(String outFileName) throws Exception{
        tstSnapshot(outFileName, true);
    }

    public void tstSnapshot(String outFileName, boolean shouldSucceed) throws Exception{

        File out = new File(snapshotDir, outFileName);
        assertFalse(out.exists());

        exec.execute("snapshot " + outFileName);

        assertEquals(shouldSucceed, out.exists());
    }

    @Test
    public void testLoadURL() throws Exception{
        String urlPath = "ftp://ftp.broadinstitute.org/distribution/igv/TEST/cpgIslands%20with%20spaces.hg18.bed";
        exec.loadFiles(urlPath, null, true, "hasSpaces");

        String localPath = TestUtils.DATA_DIR + "bed/test.bed";
        exec.loadFiles(localPath, null, true, null);

        assertEquals(2, igv.getAllTracks().size());
    }

}
