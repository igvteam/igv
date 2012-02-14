/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTIES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.ui;

import org.broad.igv.session.IGVSessionReader;
import org.broad.igv.session.SessionReader;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.Track;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.*;

import java.io.BufferedInputStream;
import java.io.InputStream;
import java.util.List;

import static org.junit.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012/02/08
 */
public class IGVTest {

    private static IGV igv;

    public IGVTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        igv = TestUtils.startGUI();
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
        TestUtils.stopGUI();
        igv = null;
    }

    @Before
    public void setUp() {

    }

    @After
    public void tearDown() {
    }

    @Test
    public void testSorts() throws Exception {
        String sessionPath = TestUtils.DATA_DIR + "/sessions/BRCA_loh2.xml";

        InputStream inputStream = new BufferedInputStream(ParsingUtils.openInputStreamGZ(new ResourceLocator(sessionPath)));
        final SessionReader sessionReader = new IGVSessionReader(igv);
        sessionReader.loadSession(inputStream, igv.getSession(), sessionPath);

        RegionScoreType[] types = RegionScoreType.values();
        int count = 0;
        for (RegionScoreType type : types) {
            count += tstSort(type);

        }
        assertTrue("Did not check enough tracks", count > 2 * types.length);
    }

    private int tstSort(RegionScoreType type) throws Exception {
        //Sort the "sortable" tracks
        ReferenceFrame frame = FrameManager.getDefaultFrame();
        final int zoom = Math.max(0, frame.getZoom());
        final String chr = frame.getChrName();
        final int start = (int) frame.getOrigin();
        final int end = (int) frame.getEnd();

        igv.sortByRegionScore(null, type, frame);

        //Need to wait for GUI to repaint
        Thread.sleep(500);

        //Check sort
        List<Track> tracks = igv.getAllTracks(false);

        Track lastTrack = null;
        int count = 0;
        for (int ii = 0; ii < tracks.size(); ii++) {
            Track track = tracks.get(ii);
            if (track.isRegionScoreType(type)) {
                String name = track.getName().toLowerCase();
                if (name.contains("reference")
                        || name.contains("refseq")) {
                    Thread.sleep(100);
                    continue;
                }
                count++;
                if (lastTrack == null) {
                    lastTrack = track;
                    continue;
                }

                float s2 = track.getRegionScore(chr, start, end, zoom, type, frame);
                float s1 = lastTrack.getRegionScore(chr, start, end, zoom, type, frame);
                assertTrue("Track named " + track.getName() + ", " + s2 + " and " + lastTrack.getName() + ", " + s1 + " out of order type " + type,
                        s2 >= s1);

                lastTrack = track;
            }
        }

        return count;
    }
}
