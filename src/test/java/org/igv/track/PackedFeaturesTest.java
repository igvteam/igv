package org.igv.track;

import htsjdk.tribble.Feature;
import org.igv.AbstractHeadlessTest;
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


        PackedFeatures<TestFeature> pf = new PackedFeatures("chr1", 0, 1000, features.iterator(), Track.DisplayMode.EXPANDED,  false );
        assertEquals(4, pf.getRowCount());

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
