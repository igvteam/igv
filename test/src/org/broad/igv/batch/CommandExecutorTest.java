/*
 * Copyright (c) 2007-2012 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTIES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NON-INFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
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

    }

    @After
    public void tearDown() throws Exception {
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
                assertEquals(maxDepth, ((AlignmentTrack) track).getMaxDepth());
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


        assertEquals(md, track1.getMaxDepth());
        assertEquals(50, track2.getMaxDepth());
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
