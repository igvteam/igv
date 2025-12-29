package org.igv.variant.vcf;


import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Jim Robinson
 * @date Aug 1, 2011
 */
public class VCFGenotype implements org.igv.variant.Genotype {

    Genotype vcfGenotype;
    Map<String, Object> attributes;
    List<org.igv.variant.Allele> alleles;

    public VCFGenotype(Genotype wrappedGenotype) {
        vcfGenotype = wrappedGenotype;
        attributes = new HashMap<String, Object>();


        if (vcfGenotype.hasDP()) {
            attributes.put("DP", vcfGenotype.getDP());
        }

        if (vcfGenotype.hasAD()) {
            String adString = "";
            int[] ad = vcfGenotype.getAD();
            for (int i = 0; i < ad.length; i++) {
                if (i != 0) adString += ",";
                adString += String.valueOf(ad[i]);
            }
            attributes.put("AD", adString);
        }

        if (vcfGenotype.hasGQ()) {
            attributes.put("GQ", vcfGenotype.getGQ());
        }

        if (vcfGenotype.hasPL()) {
            String plString = "";
            int[] pl = vcfGenotype.getPL();
            for (int i = 0; i < pl.length; i++) {
                if (i != 0) plString += ",";
                plString += String.valueOf(pl[i]);

            }
            attributes.put("PL", plString);
        }


        attributes.putAll(vcfGenotype.getExtendedAttributes());
    }

    public Map<String, Object> getAttributes() {
        return attributes;
    }

    public String getGenotypeString() {
        return vcfGenotype.getGenotypeString();
    }

    public double getPhredScaledQual() {
        return vcfGenotype.getGQ();
    }

    public String getTypeString() {
        return vcfGenotype.getType().toString();
    }

    public GenotypeType getType() {
        return vcfGenotype.getType();
    }

    public boolean isHomVar() {
        return vcfGenotype.isHomVar();
    }

    public boolean isHet() {
        return vcfGenotype.isHet();
    }

    public boolean isHomRef() {
        return vcfGenotype.isHomRef();
    }

    public boolean isNoCall() {
        return vcfGenotype.isNoCall();
    }

    public List<org.igv.variant.Allele> getAlleles() {
        if (alleles == null) {
            alleles = new ArrayList<org.igv.variant.Allele>();
            List<htsjdk.variant.variantcontext.Allele> tmp = vcfGenotype.getAlleles();
            if (tmp != null) {
                for (Allele a : tmp) {
                    alleles.add(new VCFAllele(a));
                }
            }
        }
        return alleles;
    }

    public double getAttributeAsDouble(String s) {
        return vcfGenotype.getAttributeAsDouble(s, Double.NaN);
    }
}
