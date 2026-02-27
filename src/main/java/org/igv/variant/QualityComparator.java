package org.igv.variant;

import java.util.Comparator;

class QualityComparator implements Comparator<String> {

    Variant variant;
    boolean direction;

    QualityComparator(Variant variant, boolean direction) {
        this.variant = variant;
        this.direction = direction;
    }

    public int compare(String s1, String s2) {

        double qual1 = variant.getGenotype(s1).getPhredScaledQual();
        double qual2 = variant.getGenotype(s2).getPhredScaledQual();

        int sign = direction ? -1 : 1;
        return sign * Double.compare(qual1, qual2);

    }
}
