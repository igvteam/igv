/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.variant.vcf;


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
            List<htsjdk.variant.variantcontext.Allele> tmp = vcfGenotype.getAlleles();
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
