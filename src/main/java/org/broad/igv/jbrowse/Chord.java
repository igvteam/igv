package org.broad.igv.jbrowse;

import org.broad.igv.bedpe.BedPEFeature;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.ReadMate;
import org.broad.igv.variant.Variant;

import java.util.Map;

class Chord {
    String uniqueId;
    String color;
    String refName;
    int start;
    int end;
    Mate mate;

    public Chord(BedPEFeature f) {
        this.uniqueId = f.chr1 + ":" + f.start1 + "-" + f.end1 + "_" + f.chr2 + ":" + f.start2 + "-" + f.end2;
        this.refName = f.chr1.startsWith("chr") ? f.chr1.substring(3) : f.chr1;
        this.start = f.start1;
        this.end = f.end1;
        this.mate = new Mate(f.chr2.startsWith("chr") ? f.chr2.substring(3) : f.chr2, f.start2, f.end2);
        this.color = "rgba(0, 0, 255, 0.1)";
    }

    public Chord(Alignment a) {
        ReadMate mate = a.getMate();
        this.uniqueId = a.getReadName();
        this.refName = a.getChr().startsWith("chr") ? a.getChr().substring(3) : a.getChr();
        this.start = a.getStart();
        this.end = a.getEnd();
        this.mate = new Mate(mate.getChr().startsWith("chr") ? mate.getChr().substring(3) : mate.getChr(),
                mate.getStart(), mate.getStart() + 1);
        this.color = "rgba(0, 0, 255, 0.02)";
    }

    public Chord(Variant v) {

        Map<String, Object> attrs = v.getAttributes();
        String chr2 = shortName(attrs.get("CHR2").toString());
        int end2 = Integer.parseInt(attrs.get("END").toString());
        int start2 = end2 - 1;
        String chr1 = shortName(v.getChr());
        int start1 = v.getStart();
        int end1 = v.getEnd();

        this.uniqueId = chr1 + "_" + start1 + ":" + end1 + "-" + chr2 + "_" + start2 + ":" + end2;
        this.refName = chr1;
        this.start = start1;
        this.end = end1;
        this.mate = new Mate(chr2, start2, end2);
        this.color = "rgb(0,0,255)";
    }

    static String shortName(String chr) {
        return chr.startsWith("chr") ? chr.substring(3) : chr;
    }

}
