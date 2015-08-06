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

package org.broad.igv.ui;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.igv.util.Utilities;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012/02/08
 */
public class IGVTestHeadless extends AbstractHeadlessTest{


    @Rule
    public TestRule testTimeout = new Timeout((int) 60e3);

    public IGVTestHeadless() {
    }

    @Test
    public void testSorts() throws Exception {
        String sessionPath = TestUtils.DATA_DIR + "sessions/metabric_expression.xml";
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
        final int start = Integer.parseInt("14104912") - 1;
        final int end = Integer.parseInt("36031032");
        RegionOfInterest roi = new RegionOfInterest(chr, start, end, "");

        IGV.sortByRegionScore(tracks, roi, type, frame);
        return checkIsSorted(tracks, roi, type, zoom, frame);
    }

    public static int checkIsSorted(List<Track> tracks, RegionOfInterest roi, RegionScoreType type, int zoom, ReferenceFrame frame) {

        String chr = roi.getChr();
        int start = roi.getStart();
        int end = roi.getEnd();
        String frameName = frame != null ? frame.getName() : null;

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
                float s2 = track.getRegionScore(chr, start, end, zoom, type, frameName);
                float s1 = lastTrack.getRegionScore(chr, start, end, zoom, type, frameName);
                assertTrue("Track named " + track.getName() + ", " + s2 + " and " + lastTrack.getName() + ", " + s1 + " out of order type " + type, s2 <= s1);

                lastTrack = track;
            }
        }

        return count;
    }

}
