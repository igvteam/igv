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

package org.broad.igv.track;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.util.TestUtils;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author jrobinso
 * @date Jun 10, 2011
 */
public class PackedFeaturesTest extends AbstractHeadlessTest {

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


    private void assertPackedFeaturesEqual(PackedFeatures<? extends Feature> expected, PackedFeatures<? extends Feature> actual) {

        TestUtils.assertFeatureListsEqual(expected.getFeatures().iterator(), actual.getFeatures().iterator());
        assertTrue(expected.getFeatures().size() > 0);
        assertEquals(expected.getMaxFeatureLength(), actual.getMaxFeatureLength());
        assertEquals(expected.getRowCount(), actual.getRowCount());
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

        @Override
        public String getContig() {
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
