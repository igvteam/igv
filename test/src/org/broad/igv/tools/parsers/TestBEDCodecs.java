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

package org.broad.igv.tools.parsers;

import org.broad.igv.Globals;
import org.broad.igv.feature.tribble.IGVBEDCodec;
import org.broad.igv.tools.IgvTools;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.bed.BEDCodec;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.util.Iterator;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012/01/17
 */
public class TestBEDCodecs {
    IgvTools igvTools;

    @Before
    public void setUp() throws Exception {
        Globals.setHeadless(true);
        igvTools = new IgvTools();
    }

    @After
    public void tearDown() throws Exception {
        igvTools = null;
    }

    @Test
    public void testIntervalTest() throws Exception {
        intervalTestFile(new BEDCodec());
        intervalTestFile(new IGVBEDCodec());
    }


    @Test
    public void testLargeBedNoHeader() throws Exception {
        String bedFile = TestUtils.DATA_DIR + "/bed/Unigene.noheader.sorted.bed";
        testUnigeneBed(bedFile, new BEDCodec());

        testUnigeneBed(bedFile, new IGVBEDCodec());
    }

    @Test
    public void testLargeBedWithHeader() throws Exception {
        String bedFile = TestUtils.DATA_DIR + "/bed/Unigene.withheader.sorted.bed";
        testUnigeneBed(bedFile, new BEDCodec());

        testUnigeneBed(bedFile, new IGVBEDCodec());
    }

    @Test
    public void testLargeBedWeirdHeader() throws Exception {
        String bedFile = TestUtils.DATA_DIR + "/bed/Unigene.weirdheader.sorted.bed";
        testUnigeneBed(bedFile, new BEDCodec());

        testUnigeneBed(bedFile, new IGVBEDCodec());
    }

    @Test
    public void testLargeBedNoTrack() throws Exception {
        String bedFile = TestUtils.DATA_DIR + "/bed/Unigene.notrack.sorted.bed";
        testUnigeneBed(bedFile, new BEDCodec());

        testUnigeneBed(bedFile, new IGVBEDCodec());
    }

    public void intervalTestFile(FeatureCodec codec) throws Exception {
        int startOffset = 0;
        if (codec instanceof BEDCodec) {
            startOffset = ((BEDCodec) codec).getStartOffset();
        }

        String bedFile = TestUtils.DATA_DIR + "/bed/intervalTest.bed";

        // Interval index
        TestUtils.createIndex(bedFile, IgvTools.INTERVAL_INDEX, 100);

        AbstractFeatureReader bfr = AbstractFeatureReader.getFeatureReader(bedFile, codec);
        Iterator<Feature> iter = bfr.query("chr1", 0, Integer.MAX_VALUE);
        int count = 0;
        while (iter.hasNext()) {
            Feature feat = iter.next();
            if (count == 0) {
                //Check we don't skip first line
                assertEquals(1677535, feat.getStart() - startOffset);
            }
            count++;
        }
        assertEquals(6, count);
    }


    public void testUnigeneBed(String bedFile, FeatureCodec codec) throws Exception {
        //chr2:178,599,764-179,830,787 <- CONTAINS TTN
        int startOffset = 0;
        if (codec instanceof BEDCodec) {
            startOffset = ((BEDCodec) codec).getStartOffset();
        }

        // Interval index
        TestUtils.createIndex(bedFile, IgvTools.LINEAR_INDEX, 10000);


        String chr = "chr2";
        int start = 178707289 / 2;
        int end = 179973464 * 2;

        AbstractFeatureReader bfr = AbstractFeatureReader.getFeatureReader(bedFile, codec);
        Iterator<Feature> iter = bfr.query(chr, start, end);
        int count = 0;
        while (iter.hasNext()) {
            Feature feat = iter.next();
            if (count == 0) {
                //Check we don't skip first line
                assertEquals(178707289, feat.getStart() - startOffset);
            }
            check_feat_unigene(feat, chr, start, end);
            count++;
        }

        assertEquals(71, count);

        //Re-query with some restrictions
        count = 0;
        start = 178709699;
        end = 179721089;
        iter = bfr.query(chr, start, end);
        while (iter.hasNext()) {
            Feature feat = iter.next();
            check_feat_unigene(feat, chr, start, end);
            count++;
        }

        assertEquals(65, count);
    }

    /**
     * Some checking of features queried from a given test file
     *
     * @param feat
     * @param chr
     * @param start
     * @param end
     */
    private void check_feat_unigene(Feature feat, String chr, int start, int end) {
        assertEquals(chr, feat.getChr());
        assertTrue(feat.getEnd() > feat.getStart());
        assertTrue("Start out of range", feat.getStart() >= start);
        assertTrue("end out of range", feat.getStart() <= end);
    }
}
