package org.broad.igv.sam;

import htsjdk.samtools.TextCigarCodec;
import org.junit.Assert;
import org.junit.Test;

import java.util.List;

public class ClippingCountsTest {

    @Test
    public void testFactories() {
        List<Boolean> booleans = List.of(true, false);
        for (boolean lh : booleans) {
            for (boolean ls : booleans) {
                for (boolean rs : booleans) {
                    for (boolean rh : booleans) {
                        StringBuilder cigar = new StringBuilder();
                        cigar.append(lh ? "1H" : "");
                        cigar.append(ls ? "2S" : "");
                        cigar.append("5M");
                        cigar.append(rs ? "3S" : "");
                        cigar.append(rh ? "4H" : "");
                        ClippingCounts expected = new ClippingCounts(lh ? 1 : 0, ls ? 2 : 0, rs ? 3 : 0, rh ? 4 : 0);
                        Assert.assertEquals("from string", expected, ClippingCounts.fromCigarString(cigar.toString()));
                        Assert.assertEquals("from cigar", expected, ClippingCounts.fromCigar(TextCigarCodec.decode(cigar.toString())));
                    }
                }
            }
        }
    }

    @Test
    public void testForWiringIssues(){
        final ClippingCounts counts = new ClippingCounts(1, 2, 3, 4);
        Assert.assertEquals(1, counts.getLeftHard());
        Assert.assertEquals(2, counts.getLeftSoft());
        Assert.assertEquals(3, counts.getRightSoft());
        Assert.assertEquals(4, counts.getRightHard());
        Assert.assertTrue(counts.isClipped());
        Assert.assertTrue(counts.isLeftClipped());
        Assert.assertTrue(counts.isRightClipped());
        Assert.assertEquals(3, counts.getLeft());
        Assert.assertEquals(7, counts.getRight());

        final ClippingCounts empty = new ClippingCounts(0, 0, 0 ,0);
        Assert.assertFalse(empty.isClipped());
        Assert.assertFalse(empty.isRightClipped());
        Assert.assertFalse(empty.isLeftClipped());
        Assert.assertEquals(0,empty.getLeft());
        Assert.assertEquals(0,empty.getRight());
    }
}