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

        public int getStart() {
            return start;
        }

        public int getEnd() {
            return end;
        }
    }

}
