package org.broad.igv.bedpe;

public class HicFeature extends BedPEFeature {

    private float counts;
    private float value;

    public HicFeature(String chr1, int start1, int end1, String chr2, int start2, int end2, float counts, float value) {
        super(chr1, start1, end1, chr2, start2, end2);
        this.counts = counts;
        this.value = value;
    }

    public float getValue() {
        return value;
    }

    @Override
    public String getValueString() {
        StringBuffer buf = new StringBuffer(super.getValueString());
        buf.append("<hr>");
        buf.append(String.format("counts=%.0f", counts));
        buf.append("<br>");
        buf.append(String.format("value=%.2f", value));
        return buf.toString();
    }
}
