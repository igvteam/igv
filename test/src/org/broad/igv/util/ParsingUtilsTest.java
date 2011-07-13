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

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

import org.junit.Test;

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
    public void testParseTrackLine
            () {
        // Add your code here
    }

    @Test
    public void testGetContentLengthFTP() {
        String url = "ftp://ftp.broadinstitute.org/pub/igv/genomes/genomes.txt";
        assertTrue(ParsingUtils.getContentLength(url) > 0);
    }
}

