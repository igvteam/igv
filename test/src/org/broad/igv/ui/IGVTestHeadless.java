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

package org.broad.igv.ui;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.igv.util.Utilities;
import org.fest.swing.fixture.ComponentFixture;
import org.fest.swing.fixture.FrameFixture;
import org.fest.swing.fixture.JComboBoxFixture;
import org.fest.swing.fixture.JPanelFixture;
import org.junit.*;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012/02/08
 */
public class IGVTestHeadless extends AbstractHeadlessTest{

    public IGVTestHeadless() {
    }

    @Test
    public void testSorts() throws Exception {
        String sessionPath = TestUtils.DATA_DIR + "sessions/BRCA_loh2.xml";
        List<Track> tracks = new ArrayList<Track>();
        InputStream cbioStream = ParsingUtils.openInputStream(sessionPath);
        Document document = Utilities.createDOMDocumentFromXmlStream(cbioStream);
        NodeList nodeTracks = document.getElementsByTagName("Resource");

        TrackLoader loader = new TrackLoader();

        for (int nt = 0; nt < nodeTracks.getLength(); nt++) {
            Node node = nodeTracks.item(nt);
            ResourceLocator locator = new ResourceLocator(node.getAttributes().getNamedItem("path").getTextContent());
            tracks.addAll(loader.load(locator, genome));
        }

        RegionScoreType[] types = RegionScoreType.values();
        int count = 0;
        for (RegionScoreType type : types) {
            count += tstSort(tracks, type);

        }
        assertTrue("Did not check enough tracks", count > 2 * types.length);
    }

    private int tstSort(List<Track> tracks, RegionScoreType type) throws Exception {
        //Sort the "sortable" tracks
        ReferenceFrame frame = null;
        final int zoom = 0;
        final String chr = "chr20";
        final int start = Integer.parseInt("14104912");
        final int end = Integer.parseInt("36031032");
        RegionOfInterest roi = new RegionOfInterest(chr, start, end, "");

        IGV.sortByRegionScore(tracks, roi, type, frame);
        return checkIsSorted(tracks, roi, type, zoom);
    }

    public static int checkIsSorted(List<Track> tracks, RegionOfInterest roi, RegionScoreType type, int zoom) {

        String chr = roi.getChr();
        int start = roi.getStart();
        int end = roi.getEnd();

        Track lastTrack = null;
        int count = 0;
        for (int ii = 0; ii < tracks.size(); ii++) {
            Track track = tracks.get(ii);
            if (track.isRegionScoreType(type)) {
                String name = track.getName().toLowerCase();
                if (name.contains("reference")
                        || name.contains("refseq")) {
                    continue;
                }
                count++;
                if (lastTrack == null) {
                    lastTrack = track;
                    continue;
                }

                // Test sort order -- by default tracks should be sorted in descending value
                float s2 = track.getRegionScore(chr, start, end, zoom, type, null);
                float s1 = lastTrack.getRegionScore(chr, start, end, zoom, type, null);
                assertTrue("Track named " + track.getName() + ", " + s2 + " and " + lastTrack.getName() + ", " + s1 +
                        " out of order type " + type, s2 <= s1);

                lastTrack = track;
            }
        }

        return count;
    }

}
