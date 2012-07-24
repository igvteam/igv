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


import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * @author Jim Robinson
 * @date Aug 1, 2011
 */
public class VCFGenotype implements org.broad.igv.variant.Genotype {

    Genotype vcfGenotype;
    List<org.broad.igv.variant.Allele> alleles;

    public VCFGenotype(Genotype wrappedGenotype) {
        vcfGenotype = wrappedGenotype;
    }

    public String getAttributeAsString(String key) {
        Object x = vcfGenotype.getExtendedAttribute(key);
        if (x == null) return null;
        if (x instanceof String) return (String) x;
        return String.valueOf(x); // throws an exception if this isn't a string
    }

    public Map<String, Object> getAttributes() {
        return vcfGenotype.getExtendedAttributes();
    }

    public String getGenotypeString() {
        return vcfGenotype.getGenotypeString();
    }

    public double getPhredScaledQual() {
        return vcfGenotype.getGQ();
    }

    public String getType() {
        return vcfGenotype.getType().toString();
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
            List<Allele> tmp = vcfGenotype.getAlleles();
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
