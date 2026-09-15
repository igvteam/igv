package org.igv.alignment.fiberseq;

import org.igv.alignment.AlignmentBlock;
import org.igv.alignment.AlignmentBlockImpl;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.*;

public class FiberseqAnnotationsTest {

    private static final byte[] BASES = "ACGTACGTACGTACGTACGT".getBytes();   // read length 20

    private static AlignmentBlock block(int refStart, int readOffset, int length) {
        return new AlignmentBlockImpl(refStart, BASES, null, readOffset, length, 'M');
    }

    private static AlignmentBlock softClip(int refStart, int readOffset, int length) {
        AlignmentBlockImpl b = new AlignmentBlockImpl(refStart, BASES, null, readOffset, length, 'S');
        b.setSoftClipped(true);
        return b;
    }

    @Test
    public void forwardStrand() {
        AlignmentBlock[] blocks = {block(100, 0, 20)};
        FiberseqAnnotations a = FiberseqAnnotations.create(new int[]{2}, new int[]{5}, null, null, null, 20, false, blocks);
        assertEquals(List.of(new FiberseqAnnotations.Interval(102, 107, 0)), a.getNucleosomes());
        assertTrue(a.getMsps().isEmpty());
    }

    @Test
    public void reverseStrandIsFlippedAcrossReadLength() {
        // Molecular [2,7) is [13,18) of the stored sequence for a reverse-strand read of length 20
        AlignmentBlock[] blocks = {block(100, 0, 20)};
        FiberseqAnnotations a = FiberseqAnnotations.create(new int[]{2}, new int[]{5}, null, null, null, 20, true, blocks);
        assertEquals(List.of(new FiberseqAnnotations.Interval(113, 118, 0)), a.getNucleosomes());
    }

    @Test
    public void reverseStrandWithLeadingSoftClip() {
        // 3S17M aligned at 100: read index 3 is reference 100, so read [13,18) is reference [110,115)
        AlignmentBlock[] hidden = {block(100, 3, 17)};
        AlignmentBlock[] shown = {softClip(97, 0, 3), block(100, 3, 17)};
        for (AlignmentBlock[] blocks : List.of(hidden, shown)) {
            FiberseqAnnotations a = FiberseqAnnotations.create(new int[]{2}, new int[]{5}, null, null, null, 20, true, blocks);
            assertEquals(List.of(new FiberseqAnnotations.Interval(110, 115, 0)), a.getNucleosomes());
        }
    }

    @Test
    public void intervalSpansDeletion() {
        // 10M5D10M at 100: read [8,12) covers reference 108-109 and 115-116
        AlignmentBlock[] blocks = {block(100, 0, 10), block(115, 10, 10)};
        FiberseqAnnotations a = FiberseqAnnotations.create(new int[]{8}, new int[]{4}, null, null, null, 20, false, blocks);
        assertEquals(List.of(new FiberseqAnnotations.Interval(108, 117, 0)), a.getNucleosomes());
    }

    @Test
    public void mspQualityIsUnsigned() {
        AlignmentBlock[] blocks = {block(100, 0, 20)};
        FiberseqAnnotations a = FiberseqAnnotations.create(null, null,
                new int[]{0, 10}, new int[]{4, 5}, new byte[]{0, (byte) 200}, 20, false, blocks);
        assertEquals(List.of(new FiberseqAnnotations.Interval(100, 104, 0), new FiberseqAnnotations.Interval(110, 115, 200)),
                a.getMsps());
    }

    @Test
    public void unusableIntervalsAreSkipped() {
        // Interval entirely within a soft clip, one past the read end, and mismatched MSP arrays
        AlignmentBlock[] blocks = {softClip(97, 0, 3), block(100, 3, 17)};
        FiberseqAnnotations a = FiberseqAnnotations.create(new int[]{0, 18, 5}, new int[]{3, 5, 2},
                new int[]{1, 2}, new int[]{1}, null, 20, false, blocks);
        assertEquals(List.of(new FiberseqAnnotations.Interval(102, 104, 0)), a.getNucleosomes());
        assertTrue(a.getMsps().isEmpty());
    }

    @Test
    public void noTagsReturnsNull() {
        AlignmentBlock[] blocks = {block(100, 0, 20)};
        assertNull(FiberseqAnnotations.create(null, null, null, null, null, 20, false, blocks));
    }
}
