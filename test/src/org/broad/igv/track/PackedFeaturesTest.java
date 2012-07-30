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

package org.broad.igv.track;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author jrobinso
 * @date Jun 10, 2011
 */
public class PackedFeaturesTest extends AbstractHeadlessTest{

    /**
     * Simple test of packing a list of 4 overlapping features.  The row count after packing should be 4.
     *
     * @throws Exception
     */
    @Test
    public void testGetRows() throws Exception {

        List<TestFeature> features = Arrays.asList(
                new TestFeature("chr1", 1, 100),
                new TestFeature("chr1", 2, 101),
                new TestFeature("chr1", 3, 102),
                new TestFeature("chr1", 4, 103));


        PackedFeatures<TestFeature> pf = new PackedFeatures("chr1", 0, 1000, features.iterator(), "");
        assertEquals(4, pf.getRowCount());

    }

    @Test
    public void testMerge() throws Exception{

        String filePath = TestUtils.DATA_DIR + "bed/test.bed";
        TestUtils.createIndex(filePath);
        String chr = "chr1";
        int start1 = 0;
        int end1 = 250;

        int start2 = 199;
        int end2 = 100010;

        FeatureCodec codec = CodecFactory.getCodec(filePath, genome);
        AbstractFeatureReader<Feature> bfs = AbstractFeatureReader.getFeatureReader(filePath, codec, true);
        Iterable<Feature> iter1 = bfs.query(chr, start1, end1);
        Iterable<Feature> iter2 = bfs.query(chr, start2, end2);
        Iterable<Feature> iterExp = bfs.query(chr, start1, end2);

        PackedFeatures<Feature> pFeatsMerged = new PackedFeatures<Feature>(chr, start1, end1, iter1.iterator(), "iter1");
        PackedFeatures<Feature> pFeats2 = new PackedFeatures<Feature>(chr, start2, end2, iter2.iterator(), "iter2");
        PackedFeatures<Feature> pFeatsExp = new PackedFeatures<Feature>(chr, start1, end2, iterExp.iterator(), "iterExp");

        pFeatsMerged.merge(pFeats2);

        assertPackedFeaturesEqual(pFeatsExp, pFeatsMerged);
    }

    private void assertPackedFeaturesEqual(PackedFeatures<? extends Feature> expected, PackedFeatures<? extends Feature> actual){

        TestUtils.assertFeatureListsEqual(expected.getFeatures().iterator(), actual.getFeatures().iterator());
        assertTrue(expected.getFeatures().size() > 0);
        assertEquals(expected.getMaxFeatureLength(), actual.getMaxFeatureLength());
        assertEquals(expected.getRowCount(), actual.getRowCount());
    }

    @Test
    public void testTrimTo() throws Exception{

        String filePath = TestUtils.DATA_DIR + "bed/test.bed";
        TestUtils.createIndex(filePath);
        String chr = "chr1";
        int start1 = 0;
        int end1 = 100010;

        int start2 = 199;
        int end2 = 250;

        FeatureCodec codec = CodecFactory.getCodec(filePath, genome);
        AbstractFeatureReader<Feature> bfs = AbstractFeatureReader.getFeatureReader(filePath, codec, true);
        Iterable<Feature> iter1 = bfs.query(chr, start1, end1);
        Iterable<Feature> iter2 = bfs.query(chr, start2, end2);

        PackedFeatures<Feature> pFeatsTrimmed = new PackedFeatures<Feature>(chr, start1, end1, iter1.iterator(), "iter1");
        PackedFeatures<Feature> pFeatsExp = new PackedFeatures<Feature>(chr, start2, end2, iter2.iterator(), "iter2");

        pFeatsTrimmed.trimTo(chr, start2, end2, -1);

        assertPackedFeaturesEqual(pFeatsExp, pFeatsTrimmed);
    }

    static class TestFeature implements Feature {
        String chr;
        int start;
        int end;

        TestFeature(String chr, int start, int end) {
            this.chr = chr;
            this.start = start;
            this.end = end;
        }

        public String getChr() {
            return chr;
        }

        public int getStart() {
            return start;
        }

        public int getEnd() {
            return end;
        }
    }

}
