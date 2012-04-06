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

import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.IGVTest;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.*;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import static junit.framework.Assert.*;

/**
 * User: jacob
 * Date: 2012/03/21
 */
public class CommandExecutorTest {
    static IGV igv;
    CommandExecutor exec = new CommandExecutor();

    @BeforeClass
    public static void setUpClass() throws Exception {
        igv = TestUtils.startGUI();
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
        TestUtils.stopGUI();
    }

    @Before
    public void setUp() throws Exception {
        igv.getSession().clearRegionsOfInterest();
    }

    @After
    public void tearDown() throws Exception {
        igv.removeTracks(igv.getAllTracks(false));
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
        String file = TestUtils.LARGE_DATA_DIR + "/HG00171.hg18.bam";
        igv.loadResources(Arrays.asList(new ResourceLocator(file)));

        setCheckMaxDepth(5);
        setCheckMaxDepth(50);
    }

    private void setCheckMaxDepth(int maxDepth) {
        String res = exec.execute("maxDepth " + maxDepth);
        assertFalse(res.contains("UNKNOWN"));

        List<Track> curTracks = igv.getAllTracks(false);
        int checked = 0;
        for (Track track : curTracks) {
            if (track instanceof AlignmentTrack) {
                assertEquals(maxDepth, ((AlignmentTrack) track).getDownsampleCount());
                checked++;
            }
        }
        assertTrue("No alignmenttracks found", checked > 0);
    }

    @Test
    public void testSetMaxDepthByTrack() throws Exception {

        String name1 = "HG00171.hg18.bam";
        String name2 = "HG00171.hg18.sam";
        String file1 = TestUtils.LARGE_DATA_DIR + File.separator + name1;
        String file2 = TestUtils.LARGE_DATA_DIR + File.separator + name2;
        List<ResourceLocator> resources = Arrays.asList(new ResourceLocator(file1), new ResourceLocator(file2));
        igv.loadResources(resources);

        setCheckMaxDepth(5);
        setCheckMaxDepth(50);

        int md = 100;
        exec.execute("maxDepth " + md + " " + name1);

        List<AlignmentTrack> tracks = getAlTracks(igv.getAllTracks(false));
        AlignmentTrack track1 = tracks.get(0);
        AlignmentTrack track2 = tracks.get(1);


        assertEquals(md, track1.getDownsampleCount());
        assertEquals(50, track2.getDownsampleCount());
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
        String sessionPath = TestUtils.DATA_DIR + "/sessions/BRCA_loh2.xml";
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

                tracks = igv.getAllTracks(false);
                IGVTest.checkIsSorted(tracks, roi, type, FrameManager.getDefaultFrame().getZoom());
                count++;
            }
        }


    }

}
