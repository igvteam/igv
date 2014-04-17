/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.variant.vcf;


import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypeType;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Jim Robinson
 * @date Aug 1, 2011
 */
public class VCFGenotype implements org.broad.igv.variant.Genotype {

    Genotype vcfGenotype;
    Map<String, Object> attributes;
    List<org.broad.igv.variant.Allele> alleles;

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

    public List<org.broad.igv.variant.Allele> getAlleles() {
        if (alleles == null) {
            alleles = new ArrayList<org.broad.igv.variant.Allele>();
            List<org.broadinstitute.variant.variantcontext.Allele> tmp = vcfGenotype.getAlleles();
            if (tmp != null) {
                for (Allele a : tmp) {
                    alleles.add(new VCFAllele(a.getBases()));
                }
            }
        }
        return alleles;
    }

    public double getAttributeAsDouble(String s) {
        return vcfGenotype.getAttributeAsDouble(s, Double.NaN);
    }
}
