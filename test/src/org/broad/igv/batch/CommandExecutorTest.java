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
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.*;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import static junit.framework.Assert.assertEquals;

/**
 * User: jacob
 * Date: 2012/03/21
 */
public class CommandExecutorTest {
    static IGV igv;
    static String gbmSession = TestUtils.DATA_DIR + "/sessions/gbm_subtypes.xml";
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

//        String sortType = "readGroup";
//        exec.execute("sort " + sortType);
//
//
//        List<Track> tracks = igv.getAllTracks(false);
//        System.out.println(tracks.size());
//        for(Track track: tracks){
//            if(track instanceof AlignmentTrack){
//                ((AlignmentTrack) track).sortRows();
//            }
//        }
    }

    @Test
    public void setMaxDepth() throws Exception {
        String file = TestUtils.LARGE_DATA_DIR + "/HG00171.hg18.bam";
        List<Track> loadedTracks = igv.load(new ResourceLocator(file));

        setCheckMaxDepth(5);
        setCheckMaxDepth(50);

    }

    private void setCheckMaxDepth(int maxDepth) {
        exec.execute("setMaxDepth " + maxDepth);

        List<Track> curTracks = igv.getAllTracks(false);
        for (Track track : curTracks) {
            if (track instanceof AlignmentTrack) {
                assertEquals(maxDepth, ((AlignmentTrack) track).getMaxDepth());
            }
        }
    }


}
