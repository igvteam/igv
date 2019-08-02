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
    String scoreString = "";
    double score;
    Color color;
    int thickness = 1;
    String type;
    Map<String, String> attributes;
    int row;
    BedPEShape shape;

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

    public double getMidStart() {
        return Math.min ((start1 + end1) / 2.0, (start2 + end2) / 2.0);
    }

    public double getMidEnd() {
        return Math.max ((start1 + end1) / 2.0, (start2 + end2) / 2.0);
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

    @Override
    public Color getColor() {
        return color;
    }

    @Override
    public int getThickness() {
        return thickness;
    }

    public String getContig() {
        return getChr();
    }

    @Override
    public BedPEShape getShape() {
        return shape;
    }

    @Override
    public void setShape(BedPEShape shape) {
        this.shape = shape;
    }

    public String getValueString() {

        StringBuffer buf = new StringBuffer();

        String locus1 = chr1 + ":" + start1 + "-" + end1;
        String locus2 = chr2 + ":" + start2 + "-" + end2;
        if(name != null && name.length() > 0 && !name.equals(".")) {
            buf.append(name + "<br>");
        }
        buf.append(locus1);
        buf.append("<br>" + locus2);
        buf.append("<br>Score: " + scoreString);
        if(attributes != null) {
            buf.append("<br><hr>");
            for (Map.Entry<String, String> entry : attributes.entrySet()) {
                buf.append("<br>" + entry.getKey() + ": " + entry.getValue());
            }
        }

        return buf.toString();
    }

    @Override
    public double getCenterDistance() {
        return Math.abs((start1 + end1) / 2.0 - (start2 + end2)  / 2.0);
    }
}
