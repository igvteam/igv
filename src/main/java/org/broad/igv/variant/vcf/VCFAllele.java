package org.broad.igv.variant.vcf;

import org.broad.igv.variant.Allele;

/**
 * @author Jim Robinson
 * @date Aug 1, 2011
 */
public class VCFAllele implements Allele {

    htsjdk.variant.variantcontext.Allele htsjdkAllele;


    public VCFAllele(htsjdk.variant.variantcontext.Allele htsjdkAllele) {
        this.htsjdkAllele = htsjdkAllele;
    }

    public byte[] getBases() {
        return htsjdkAllele.getBases();
    }

    public boolean isNonRefAllele() {
        return htsjdkAllele.isNonRefAllele();
    }

    @Override
    public String getDisplayString() {
        return toString();
    }

    public String toString() {
        if (htsjdkAllele.isNonRefAllele()) {
            return "NON_REF";
        } else
            return htsjdkAllele.getDisplayString();
    }

}
