package org.igv.variant;

import java.util.Comparator;

class GenotypeComparator implements Comparator<String> {

    Variant variant;
    boolean direction;

    GenotypeComparator(Variant variant, boolean direction) {
        this.variant = variant;
        this.direction = direction;
    }

    public int compare(String e1, String e2) {

        int genotype1 = classifyGenotype(variant.getGenotype(e1));
        int genotype2 = classifyGenotype(variant.getGenotype(e2));

        if (genotype2 == genotype1) {
            return 0;
        } else if (genotype2 > genotype1) {
            return direction ? 1 : -1;
        } else {
            return direction ? -1 : 1;
        }
    }


    private int classifyGenotype(Genotype genotype) {

        if (genotype.isNoCall()) {
            return direction ? 1 : 10;
        } else if (genotype.isHomVar()) {
            return 4;
        } else if (genotype.isHet()) {
            return 3;
        } else if (genotype.isHomRef()) {
            return direction ? 2 : 9;
        }
        return -1; //Unknown
    }
}
