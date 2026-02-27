package org.igv.variant;

import java.util.Comparator;

class DepthComparator implements Comparator<String> {

    Variant variant;
    boolean direction;

    DepthComparator(Variant variant, boolean direction) {
        this.variant = variant;
        this.direction = direction;
    }

    public int compare(String s1, String s2) {

        double readDepth1 = variant.getGenotype(s1).getAttributeAsDouble("DP");
        double readDepth2 = variant.getGenotype(s2).getAttributeAsDouble("DP");

        int sign = direction ? -1 : 1;
        return sign * Double.compare(readDepth1, readDepth2);

    }
}
