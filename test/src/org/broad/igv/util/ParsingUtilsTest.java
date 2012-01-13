/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.util;

import junit.framework.Assert;
import org.broad.igv.Globals;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.WindowFunction;
import org.junit.Test;

import java.awt.*;

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

/**
 * User: jrobinso
 * Date: Feb 8, 2010
 */
public class ParsingUtilsTest {
    @Test
    public void testEstimateLineCount() {
        // Add your code here
    }

    @Test
    public void testSplit1() {
        String[] tokens = new String[10];
        String blankColumnLine = "a\tb\t\td";
        int nTokens = ParsingUtils.split(blankColumnLine, tokens, '\t');
        assertEquals(4, nTokens);
        assertEquals("a", tokens[0]);
        assertEquals("b", tokens[1]);
        assertEquals("", tokens[2]);
        assertEquals("d", tokens[3]);
    }

    @Test
    public void testSplit2() {
        String[] tokens = new String[10];
        String blankColumnLine = "a\tb\t\td\t";
        int nTokens = ParsingUtils.split(blankColumnLine, tokens, '\t');
        assertEquals(5, nTokens);
        assertEquals("a", tokens[0]);
        assertEquals("b", tokens[1]);
        assertEquals("", tokens[2]);
        assertEquals("d", tokens[3]);
        assertEquals("", tokens[2]);
    }

    @Test
    public void testSplitWhitespace1() {
        String[] tokens = new String[10];
        String blankColumnLine = "a b\t\td";
        int nTokens = ParsingUtils.splitWhitespace(blankColumnLine, tokens);
        assertEquals(4, nTokens);
        assertEquals("a", tokens[0]);
        assertEquals("b", tokens[1]);
        assertEquals("", tokens[2]);
        assertEquals("d", tokens[3]);
    }

    @Test
    public void testSplitWhitespace2() {
        String[] tokens = new String[10];
        String blankColumnLine = "a b\t\td\t";
        int nTokens = ParsingUtils.splitWhitespace(blankColumnLine, tokens);
        assertEquals(5, nTokens);
        assertEquals("a", tokens[0]);
        assertEquals("b", tokens[1]);
        assertEquals("", tokens[2]);
        assertEquals("d", tokens[3]);
        assertEquals("", tokens[2]);
    }

    @Test
    public void testSplitWhitespace3() {
        String[] tokens = new String[10];
        String blankColumnLine = "a   b  \t \td";
        int nTokens = ParsingUtils.splitWhitespace(blankColumnLine, tokens);
        assertEquals(5, nTokens);
        assertEquals("a", tokens[0]);
        assertEquals("b", tokens[1]);
        assertEquals("", tokens[2]);
        assertEquals("", tokens[2]);
        assertEquals("d", tokens[4]);
    }

    @Test
    public void testComputeReadingShifts
            () {
        // Add your code here
    }

    @Test
    public void testGetContentLengthFTP() {
        assertTrue(ParsingUtils.getContentLength(TestUtils.AVAILABLE_FTP_URL) > 0);

        long start_time = System.currentTimeMillis();
        assertTrue(ParsingUtils.getContentLength(TestUtils.UNAVAILABLE_FTP_URL) == -1);
        long end_time = System.currentTimeMillis();
        assertTrue(end_time - start_time < Globals.CONNECT_TIMEOUT + 1000);
    }

    @Test
    public void testParseInt() {
        String with_commas = "123456";
        int expected = 123456;
        int actual = ParsingUtils.parseInt(with_commas);
        assertEquals(expected, actual);

        String exp_not = "3.5e4";
        expected = 35000;
        assertEquals(expected, ParsingUtils.parseInt(exp_not));
    }

    @Test
    public void testParseTrackLine() {
        String trackLine = "track type=bigWig name=\"Track 196\" visibility=2 " +
                "description=\" CD34 - H3K27me3 - hg19 - 18.7 M/20.9 M - 61P7DAAXX.6\" " +
                "maxHeightPixels=70 viewLimits=0:18 windowingFunction=mean autoScale=off " +
                "bigDataUrl=http://www.broadinstitute.org/epigenomics/dataportal/track_00196.portal.bw " +
                "color=255,0,0";

        TrackProperties props = new TrackProperties();
        ParsingUtils.parseTrackLine(trackLine, props);
        assertEquals("Track 196", props.getName());
        assertEquals(Track.DisplayMode.EXPANDED, props.getDisplayMode());
        assertEquals(" CD34 - H3K27me3 - hg19 - 18.7 M/20.9 M - 61P7DAAXX.6", props.getDescription());
        assertEquals(70, props.getHeight());
        assertEquals(0, props.getMinValue(), 1.0e-9);
        assertEquals(18, props.getMaxValue(), 1.0e-9);
        assertEquals(WindowFunction.mean, props.getWindowingFunction());
        assertEquals(false, props.isAutoScale());
        assertEquals(new Color(255, 0, 0), props.getColor());
        assertEquals("http://www.broadinstitute.org/epigenomics/dataportal/track_00196.portal.bw", props.getDataURL());
    }
}

