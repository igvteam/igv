package org.igv.circview.model;

import htsjdk.samtools.SAMTag;
import org.igv.Globals;
import org.igv.bedpe.BedPE;
import org.igv.sam.Alignment;
import org.igv.sam.ReadMate;
import org.igv.sam.SupplementaryAlignment;
import org.igv.variant.Variant;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * A chord (arc) connecting one genomic region to its {@link Mate} region.
 *
 * <p>Mirrors the chord feature objects passed to addChords() in circularView.js:
 * <pre>
 * {
 *   uniqueId, color, refName, start, end,
 *   mate: { refName, start, end }
 * }
 * </pre>
 * The per-feature {@code color} is optional; when null the owning
 * {@link ChordSet}'s color is used.
 */
public final class Chord {

    private  String uniqueId;
    private  String refName;
    private  long start;
    private  long end;
    private  Mate mate;
    private  Color color;

    public Chord() {
    }

    public Chord(String uniqueId, String refName, long start, long end, Mate mate, Color color) {
        this.uniqueId = uniqueId;
        this.refName = refName;
        this.start = start;
        this.end = end;
        this.mate = mate;
        this.color = color;
    }

    public static Chord fromBedPE(BedPE f) {

        Chord c = new Chord();
        String chr1 = shortName(f.getChr1());
        String chr2 = shortName(f.getChr2());
        c.uniqueId = chr1 + ":" + f.getStart1() + "-" + f.getEnd1() + "_" + chr2 + ":" + f.getStart2() + "-" + f.getEnd2();
        c.refName = chr1;
        c.start = f.getStart1();
        c.end = f.getEnd1();
        c.mate = new Mate(chr2, f.getStart2(), f.getEnd2());
        return c;
    }

    public static Chord fromPEAlignment(Alignment a) {
        Chord c = new Chord();
        ReadMate readMate = a.getMate();
        c.uniqueId = a.getReadName();
        c.refName = shortName(a.getChr());
        c.start = a.getStart();
        c.end = a.getEnd();
        c.mate = new Mate(shortName(readMate.getChr()),
                readMate.getStart(),
                readMate.getStart() + 1);
        return c;
    }

    public static List<Chord> fromSAString(Alignment a) {
        String sastring = a.getAttribute(SAMTag.SA.name()).toString();
        String[] records = Globals.semicolonPattern.split(sastring);
        List<Chord> chords = new ArrayList<>();
        int n = 0;
        for (String rec : records) {
            SupplementaryAlignment sa = SupplementaryAlignment.fromSingleSaTagRecord(rec);
            if (sa.start != a.getStart()) {
                Chord c = new Chord();
                c.uniqueId = a.getReadName() + "_" + n++;
                c.refName = shortName(a.getChr());
                c.start = a.getStart();
                c.end = a.getEnd();
                c.mate = new Mate(shortName(sa.chr), sa.start, sa.start + sa.getLengthOnReference());
                chords.add(c);
            }
        }
        return chords;
    }

    public static Chord fromVariant(Variant v) {
        Chord c = new Chord();
        Map<String, Object> attrs = v.getAttributes();
        String chr2 = shortName(attrs.get("CHR2").toString());
        int end2 = Integer.parseInt(attrs.get("END").toString());
        int start2 = end2 - 1;
        String chr1 = shortName(v.getChr());
        int start1 = v.getStart();
        int end1 = v.getEnd();
        c.uniqueId = chr1 + "_" + start1 + ":" + end1 + "-" + chr2 + "_" + start2 + ":" + end2;
        c.refName = chr1;
        c.start = start1;
        c.end = end1;
        c.mate = new Mate(chr2, start2, end2);
        return c;
    }

    static String shortName(String chr) {
        return chr.startsWith("chr") ? chr.substring(3) : chr;
    }
    
    public String getUniqueId() {
        return uniqueId;
    }

    public String getRefName() {
        return refName;
    }

    public long getStart() {
        return start;
    }

    public long getEnd() {
        return end;
    }

    public Mate getMate() {
        return mate;
    }

    /** Optional per-feature color; may be null (fall back to the chord set color). */
    public Color getColor() {
        return color;
    }

    @Override
    public String toString() {
        return "ChordFeature{" + refName + ":" + start + "-" + end
                + " -> " + mate.getRefName() + ":" + mate.getStart() + "-" + mate.getEnd()
                + (uniqueId != null ? " id=" + uniqueId : "") + "}";
    }
}
