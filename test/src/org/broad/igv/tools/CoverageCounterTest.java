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

package org.broad.igv.tools;

import org.broad.igv.Globals;
import org.broad.igv.tools.parsers.DataConsumer;
import org.broad.igv.track.TrackType;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.assertEquals;

/**
 * public CoverageCounter(String alignmentFile,
 * DataConsumer consumer,
 * int windowSize,
 * int extFactor,
 * File tdfFile,  // For reference
 * File wigFile,
 * Genome genome,
 * String options)
 */
public class CoverageCounterTest {

    @BeforeClass
    public static void init() {
        Globals.setHeadless(true);
    }

    @Test

    public void testMappingQualityFlag() throws IOException {
        String bamURL = "http://www.broadinstitute.org/igvdata/1KG/pilot2Bams/NA12878.SLX.bam";
        String options = "m=30,q@1:16731624-16731624";
        int windowSize = 1;

        TestDataConsumer dc = new TestDataConsumer();

        CoverageCounter cc = new CoverageCounter(bamURL, dc, windowSize, 0, null, null, null, options);

        cc.parse();

        String totalCount = dc.attributes.get("totalCount");

        assertEquals("19", totalCount);

    }

    @Test
    public void testInclueDuplicatesFlag() throws IOException {
        //TODO Guessing this URL is obsolete
        String bamURL = "http://www.broadinstitute.org/igvdata/BodyMap/IlluminaHiSeq2000_BodySites/Merged/HBM.adipose.bam.sorted.bam";
        String options = "d,q@chr1:153425249-153425249";
        int windowSize = 1;

        TestDataConsumer dc = new TestDataConsumer();

        CoverageCounter cc = new CoverageCounter(bamURL, dc, windowSize, 0, null, null, null, options);

        cc.parse();

        String totalCount = dc.attributes.get("totalCount");

        assertEquals("22", totalCount);

    }


    static class TestDataConsumer implements DataConsumer {

        Map<String, String> attributes = new HashMap();

        public void setType(String type) {
            //To change body of implemented methods use File | Settings | File Templates.
        }

        public void addData(String chr, int start, int end, float[] data, String name) {
            //To change body of implemented methods use File | Settings | File Templates.
        }

        public void parsingComplete() {
            //To change body of implemented methods use File | Settings | File Templates.
        }

        public void setTrackParameters(TrackType trackType, String trackLine, String[] trackNames) {
            //To change body of implemented methods use File | Settings | File Templates.
        }

        public void setSortTolerance(int tolerance) {
            //To change body of implemented methods use File | Settings | File Templates.
        }

        public void setAttribute(String key, String value) {
            attributes.put(key, value);
        }
    }
}
