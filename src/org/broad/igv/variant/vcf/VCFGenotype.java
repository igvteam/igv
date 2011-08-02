/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.variant.vcf;


import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;

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
        return vcfGenotype.getAttributeAsString(key);
    }

    public Map<String, Object> getAttributes() {
        return vcfGenotype.getAttributes();
    }

    public String getGenotypeString() {
        return vcfGenotype.getGenotypeString();
    }

    public double getPhredScaledQual() {
        return vcfGenotype.getPhredScaledQual();
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
            List<Allele> tmp = vcfGenotype.getAlleles();
            for (Allele a : tmp) {
                alleles.add(new VCFAllele(a.getBases()));
            }
        }
        return alleles;
    }

    public Double getAttributeAsDouble(String s) {
        return vcfGenotype.getAttributeAsDoubleNoException(s);
    }
}
