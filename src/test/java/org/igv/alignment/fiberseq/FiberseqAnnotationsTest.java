package org.igv.alignment.fiberseq;

import org.igv.alignment.AlignmentBlock;
import org.igv.alignment.AlignmentBlockImpl;
import org.igv.alignment.fiberseq.FiberseqAnnotations.Interval;
import org.junit.Test;

import java.util.List;
import java.util.Map;

import static org.junit.Assert.*;

public class FiberseqAnnotationsTest {

    private static final byte[] BASES = "ACGTACGTACGTACGTACGT".getBytes();   // read length 20

    private static final AlignmentBlock[] FULL = {block(100, 0, 20)};

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
        FiberseqAnnotations a = FiberseqAnnotations.create(new int[]{2}, new int[]{5}, null, null, null, 20, false, FULL);
        assertEquals(List.of(new Interval(102, 107, 0)), a.getNucleosomes());
        assertTrue(a.getMsps().isEmpty());
    }

    @Test
    public void reverseStrandIsFlippedAcrossReadLength() {
        // Molecular [2,7) is [13,18) of the stored sequence for a reverse-strand read of length 20
        FiberseqAnnotations a = FiberseqAnnotations.create(new int[]{2}, new int[]{5}, null, null, null, 20, true, FULL);
        assertEquals(List.of(new Interval(113, 118, 0)), a.getNucleosomes());
    }

    @Test
    public void reverseStrandWithLeadingSoftClip() {
        // 3S17M aligned at 100: read index 3 is reference 100, so read [13,18) is reference [110,115)
        AlignmentBlock[] hidden = {block(100, 3, 17)};
        AlignmentBlock[] shown = {softClip(97, 0, 3), block(100, 3, 17)};
        for (AlignmentBlock[] blocks : List.of(hidden, shown)) {
            FiberseqAnnotations a = FiberseqAnnotations.create(new int[]{2}, new int[]{5}, null, null, null, 20, true, blocks);
            assertEquals(List.of(new Interval(110, 115, 0)), a.getNucleosomes());
        }
    }

    @Test
    public void intervalSpansDeletion() {
        // 10M5D10M at 100: read [8,12) covers reference 108-109 and 115-116
        AlignmentBlock[] blocks = {block(100, 0, 10), block(115, 10, 10)};
        FiberseqAnnotations a = FiberseqAnnotations.create(new int[]{8}, new int[]{4}, null, null, null, 20, false, blocks);
        assertEquals(List.of(new Interval(108, 117, 0)), a.getNucleosomes());
    }

    @Test
    public void mspQualityIsUnsigned() {
        FiberseqAnnotations a = FiberseqAnnotations.create(null, null,
                new int[]{0, 10}, new int[]{4, 5}, new byte[]{0, (byte) 200}, 20, false, FULL);
        assertEquals(List.of(new Interval(100, 104, 0), new Interval(110, 115, 200)), a.getMsps());
    }

    @Test
    public void unusableIntervalsAreSkipped() {
        // Interval entirely within a soft clip, one past the read end, and mismatched MSP arrays
        AlignmentBlock[] blocks = {softClip(97, 0, 3), block(100, 3, 17)};
        FiberseqAnnotations a = FiberseqAnnotations.create(new int[]{0, 18, 5}, new int[]{3, 5, 2},
                new int[]{1, 2}, new int[]{1}, null, 20, false, blocks);
        assertEquals(List.of(new Interval(102, 104, 0)), a.getNucleosomes());
        assertTrue(a.getMsps().isEmpty());
    }

    @Test
    public void noTagsReturnsNull() {
        assertNull(FiberseqAnnotations.create(null, null, null, null, null, 20, false, FULL));
        assertNull(FiberseqAnnotations.fromTags(Map.of()::get, 20, false, FULL));
    }

    @Test
    public void maTagMatchesLegacyTags() {
        // Ma starts are 1-based: nucleosome 3-5 is legacy ns=2, nl=5
        Map<String, Object> ma = Map.of("Ma", "20;nuc.:3-5;msp.:11-5");
        Map<String, Object> legacy = Map.of("ns", new int[]{2}, "nl", new int[]{5}, "as", new int[]{10}, "al", new int[]{5});
        for (boolean negative : new boolean[]{false, true}) {
            FiberseqAnnotations fromMa = FiberseqAnnotations.fromTags(ma::get, 20, negative, FULL);
            FiberseqAnnotations fromLegacy = FiberseqAnnotations.fromTags(legacy::get, 20, negative, FULL);
            assertEquals(fromLegacy.getNucleosomes(), fromMa.getNucleosomes());
            assertEquals(fromLegacy.getMsps(), fromMa.getMsps());
        }
        assertEquals(List.of(new Interval(113, 118, 0)), FiberseqAnnotations.fromTags(ma::get, 20, true, FULL).getNucleosomes());
    }

    @Test
    public void maReadLengthIsTheFrameWithoutSequence() {
        FiberseqAnnotations a = FiberseqAnnotations.createFromMa("20;nuc.:3-5", null, 0, true, FULL);
        assertEquals(List.of(new Interval(113, 118, 0)), a.getNucleosomes());
    }

    @Test
    public void staleMaReadLengthIsIgnored() {
        assertNull(FiberseqAnnotations.createFromMa("21;nuc.:3-5", null, 20, false, FULL));
    }

    @Test
    public void fireQualityOverlaysMatchingMsp() {
        FiberseqAnnotations a = FiberseqAnnotations.createFromMa("20;msp.:1-4,11-5;fire.Q:11-5",
                new byte[]{(byte) 200}, 20, false, FULL);
        assertEquals(List.of(new Interval(100, 104, 0), new Interval(110, 115, 200)), a.getMsps());
    }

    @Test
    public void qualitiesAreConsumedAcrossAllSections() {
        // ctcf takes two qualities per annotation and fiberseq_callable none, so msp's quality is the fifth value.
        // A type split over several sections accumulates.
        FiberseqAnnotations a = FiberseqAnnotations.createFromMa(
                "20;ctcf+PQ:1-2,5-2;fiberseq_callable.:1-20;msp.Q:11-5;nuc-:3-2;nuc+:15-2",
                new byte[]{1, 2, 3, 4, 77}, 20, false, FULL);
        assertEquals(List.of(new Interval(110, 115, 77)), a.getMsps());
        assertEquals(List.of(new Interval(102, 104, 0), new Interval(114, 116, 0)), a.getNucleosomes());
    }

    @Test
    public void malformedMaTagIsIgnored() {
        byte[] aq = {10};
        assertNotNull(FiberseqAnnotations.createFromMa("20;msp.:11-5", aq, 20, false, FULL));
        for (String ma : List.of(
                "",
                "x;msp.:11-5",                      // read length
                "20;msp.:11-5;nuc:3-5",             // no strand
                "20;msp.:11-5;.:3-5",               // no type name
                "20;msp.:11-5;nuc.:3",              // no inline length
                "20;msp.:11-5;nuc.:0-5",            // start is 1-based
                "20;msp.:11-5;nuc.:3--5",           // negative length
                "20;msp.:11-5;nuc.X:3-5",           // unknown quality letter
                "20;msp.:11-5;fire.QQ:11-5")) {     // too few qualities
            assertNull(ma, FiberseqAnnotations.createFromMa(ma, aq, 20, false, FULL));
        }
        assertNull(FiberseqAnnotations.createFromMa("20;msp.:11-5;fire.Q:11-5", null, 20, false, FULL));
    }

    @Test
    public void tagFamilyPrecedence() {
        // Uppercase MA wins, and its qualities come from AQ, not Aq
        Map<String, Object> both = Map.of(
                "MA", "20;msp.Q:11-5", "AQ", new byte[]{50},
                "Ma", "20;msp.Q:1-4", "Aq", new byte[]{90});
        assertEquals(List.of(new Interval(110, 115, 50)), FiberseqAnnotations.fromTags(both::get, 20, false, FULL).getMsps());

        // A non-string MA does not shadow Ma
        Map<String, Object> foreignMA = Map.of("MA", 7, "Ma", "20;nuc.:3-5");
        assertEquals(List.of(new Interval(102, 107, 0)), FiberseqAnnotations.fromTags(foreignMA::get, 20, false, FULL).getNucleosomes());

        // Molecular annotation tags take precedence over legacy tags
        Map<String, Object> mixed = Map.of("Ma", "20;nuc.:3-5", "ns", new int[]{10}, "nl", new int[]{5});
        assertEquals(List.of(new Interval(102, 107, 0)), FiberseqAnnotations.fromTags(mixed::get, 20, false, FULL).getNucleosomes());

        assertTrue(FiberseqAnnotations.hasTags(foreignMA::get));
        assertTrue(FiberseqAnnotations.hasTags(Map.of("as", new int[]{1})::get));
        assertFalse(FiberseqAnnotations.hasTags(Map.of("MA", 7)::get));
    }
}
