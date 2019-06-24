package org.broad.igv.bedpe;

import java.awt.*;
import java.util.Map;

/**
 * Created by jrobinso on 6/29/18.
 */
public class BedPEFeature implements BedPE {

    String chr1;
    int start1;
    int end1;
    String chr2;
    int start2;
    int end2;
    String name;
    double score;
    Color color;
    int thickness = 1;
    String type;
    Map<String, String> attributes;
    int row;

    public BedPEFeature(String chr1, int start1, int end1, String chr2, int start2, int end2) {
        this.chr1 = chr1;
        this.start1 = start1;
        this.chr2 = chr2;
        this.end1 = end1;
        this.start2 = start2;
        this.end2 = end2;
    }

    public BedPEFeature get() {
        return this;
    }

    public String getChr() {
        if(isSameChr()) {
            return chr1;
        } else {
            return null;
        }
    }

    public int getStart() {
        return Math.min(start1, start2);
    }

    public int getEnd() {
        return Math.max(end1, end2);
    }

    @Override
    public double getScore() {
        return score;
    }

    public boolean isSameChr() {
        return chr1.equals(chr2);
    }

    @Override
    public void setRow(int row) {
        this.row = row;
    }

    @Override
    public int getRow() {
        return row;
    }

    public String getContig() {
        return getChr();
    }

}
