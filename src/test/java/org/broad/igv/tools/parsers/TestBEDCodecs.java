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

package org.broad.igv.tools.parsers;

import junit.framework.Assert;
import org.broad.igv.Globals;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.feature.tribble.FeatureFileHeader;
import org.broad.igv.feature.tribble.IGVBEDCodec;
import org.broad.igv.tools.IgvTools;
import org.broad.igv.util.TestUtils;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.bed.BEDCodec;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.util.Iterator;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertNotNull;

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
        String bedFile = TestUtils.DATA_DIR + "bed/Unigene.noheader.sorted.bed";
        tstUnigeneBed(bedFile, new BEDCodec());

        tstUnigeneBed(bedFile, new IGVBEDCodec());
    }

    @Test
    public void testLargeBedWithHeader() throws Exception {
        String bedFile = TestUtils.DATA_DIR + "bed/Unigene.withheader.sorted.bed";
        tstUnigeneBed(bedFile, new BEDCodec());

        tstUnigeneBed(bedFile, new IGVBEDCodec());
    }

    @Test
    public void testLargeBedWeirdHeader() throws Exception {
        String bedFile = TestUtils.DATA_DIR + "bed/Unigene.weirdheader.sorted.bed";
        tstUnigeneBed(bedFile, new BEDCodec());

        tstUnigeneBed(bedFile, new IGVBEDCodec());
    }

    @Test
    public void testLargeBedNoTrack() throws Exception {
        String bedFile = TestUtils.DATA_DIR + "bed/Unigene.notrack.sorted.bed";
        tstUnigeneBed(bedFile, new BEDCodec());

        tstUnigeneBed(bedFile, new IGVBEDCodec());
    }

    @Test
    public void testGffTags() throws Exception {
        String bedFile = TestUtils.DATA_DIR + "bed/gene.bed";
        Genome genome = TestUtils.loadGenome();

        FeatureCodec codec1 = CodecFactory.getCodec(bedFile, genome);
        assertTrue(codec1 instanceof IGVBEDCodec);
        IGVBEDCodec codec = (IGVBEDCodec) codec1;

        AbstractFeatureReader bfr = AbstractFeatureReader.getFeatureReader(bedFile, codec, false);
        FeatureFileHeader header = (FeatureFileHeader) bfr.getHeader();
        Assert.assertNotNull(header);
        assertTrue(codec.isGffTags());

        Iterator<BasicFeature> iter = bfr.iterator();
        while (iter.hasNext()) {
            BasicFeature feat = iter.next();
            //Note: These are not in general equal, but they are for this data file
            assertEquals(feat.getName(), feat.getIdentifier());
            assertNotNull("No ID found for feature", feat.getIdentifier());
            assertNotNull("No description found for feature", feat.getDescription());
            assertNotNull("Feature not in FeatureDB", FeatureDB.getFeature(feat.getIdentifier()));
        }


    }

    public void intervalTestFile(FeatureCodec codec) throws Exception {
        int startOffset = 0;
        if (codec instanceof BEDCodec) {
            startOffset = ((BEDCodec) codec).getStartOffset();
        }

        String bedFile = TestUtils.DATA_DIR + "bed/intervalTest.bed";

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


    public void tstUnigeneBed(String bedFile, FeatureCodec codec) throws Exception {
        //chr2:178,599,764-179,830,787 <- CONTAINS TTN
        int startOffset = 0;
        if (codec instanceof BEDCodec) {
            startOffset = ((BEDCodec) codec).getStartOffset();
        }

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


    @Test
    public void testLength1Feature() throws Exception {
        String bedFile = TestUtils.DATA_DIR + "bed/snp_calls.bed";
        TestUtils.createIndex(bedFile, IgvTools.LINEAR_INDEX, 10000);
        FeatureCodec codec = CodecFactory.getCodec(bedFile, null);

        AbstractFeatureReader<Feature, ?> bfr = AbstractFeatureReader.getFeatureReader(bedFile, codec);
        for (Feature feat : bfr.iterator()) {
            BasicFeature f = (BasicFeature) feat;
            assertEquals(1, f.getLength());
        }
    }
}
