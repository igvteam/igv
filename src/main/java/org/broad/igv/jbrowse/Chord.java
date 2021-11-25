package org.broad.igv.jbrowse;

import org.broad.igv.Globals;
import org.broad.igv.bedpe.BedPEFeature;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.ReadMate;
import org.broad.igv.sam.SupplementaryAlignment;
import org.broad.igv.variant.Variant;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

class Chord {
    String uniqueId;
    String color;
    String refName;
    int start;
    int end;
    Mate mate;

    private Chord() {
    }

    public static Chord fromBedPE(BedPEFeature f) {

        Chord c = new Chord();
        String chr1 = shortName(f.chr1);
        String chr2 = shortName(f.chr2);
        c.uniqueId = chr1 + ":" + f.start1 + "-" + f.end1 + "_" + chr2 + ":" + f.start2 + "-" + f.end2;
        c.refName = chr1;
        c.start = f.start1;
        c.end = f.end1;
        c.mate = new Mate(chr2, f.start2, f.end2);
        c.color = "rgba(0, 0, 255, 0.1)";
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
        c.color = "rgba(0, 0, 255, 0.02)";
        return c;
    }

    public static List<Chord>  fromSAString(Alignment a) {

        String sastring = a.getAttribute("SA").toString();
        String[] records = Globals.semicolonPattern.split(sastring);
        List<Chord> chords = new ArrayList<>();
        int n = 0;
        for (String rec : records) {
            SupplementaryAlignment sa = new SupplementaryAlignment(rec);
            if (sa.start != a.getStart()) {
                Chord c = new Chord();
                c.uniqueId = a.getReadName() + "_" + n++;
                c.refName = shortName(a.getChr());
                c.start = a.getStart();
                c.end = a.getEnd();
                c.mate = new Mate(shortName(sa.chr), sa.start, sa.start + sa.lenOnRef);
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
        c.color = "rgb(0,0,255)";
        return c;
    }

    public String toJson() {
        StringBuffer buf = new StringBuffer();
        buf.append("{");
        buf.append(JsonUtils.toJson("uniqueId", uniqueId));
        buf.append(",");
        buf.append(JsonUtils.toJson("color", color));
        buf.append(",");
        buf.append(JsonUtils.toJson("refName", refName));
        buf.append(",");
        buf.append(JsonUtils.toJson("start", start));
        buf.append(",");
        buf.append(JsonUtils.toJson("end", end));
        buf.append(",");
        buf.append("\"mate\":");
        buf.append(mate.toJson());
        buf.append("}");
        return buf.toString();
    }

    static String shortName(String chr) {
        return chr.startsWith("chr") ? chr.substring(3) : chr;
    }
}


class Mate {
    String refName;
    int start;
    int end;

    public Mate(String refName, int start, int end) {
        this.refName = refName;
        this.start = start;
        this.end = end;
    }

    public String toJson() {
        StringBuffer buf = new StringBuffer();
        buf.append("{");
        buf.append(JsonUtils.toJson("refName", refName));
        buf.append(",");
        buf.append(JsonUtils.toJson("start", start));
        buf.append(",");
        buf.append(JsonUtils.toJson("end", end));
        buf.append("}");
        return buf.toString();
    }
}

