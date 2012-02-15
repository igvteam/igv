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

import org.broad.igv.variant.Allele;
import org.broad.igv.variant.Genotype;
import org.broad.igv.variant.Variant;
//import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

/**
 * @author Jim Robinson
 * @date Aug 1, 2011
 */
public class VCFVariant implements Variant {

    VariantContext variantContext;
    Set<Allele> alternateAlleles;
    private ZygosityCount zygosityCount;
    private boolean isIndel;

    String chr;   // <= might override what's in file from chr alias table

    public VCFVariant(VariantContext variantContext, String chr) {
        this.variantContext = variantContext;
        this.chr = chr;
        isIndel = variantContext.getType() == VariantContext.Type.INDEL;
        init();
    }


    private void init() {
            zygosityCount = new ZygosityCount();
            for (String sample : getSampleNames()) {
                Genotype genotype = getGenotype(sample);
                zygosityCount.incrementCount(genotype);
            }
    }


    public String getID() {
        return variantContext.getID();
    }

    public boolean isFiltered() {
        return variantContext.isFiltered();
    }

    public String getAttributeAsString(String key) {
        return variantContext.getAttributeAsString(key, null);
    }

    public String getReference() {
        return variantContext.getReference().toString();
    }

    public Set<Allele> getAlternateAlleles() {
        if (alternateAlleles == null) {
            List<org.broadinstitute.sting.utils.variantcontext.Allele> tmp = variantContext.getAlternateAlleles();
            alternateAlleles = new HashSet(tmp.size());
            for (org.broadinstitute.sting.utils.variantcontext.Allele a : tmp) {
                alternateAlleles.add(new VCFAllele(a.getBases()));
            }
        }
        return alternateAlleles;
    }

    public double getPhredScaledQual() {
        return variantContext.getPhredScaledQual();
    }

    public String getType() {
        return variantContext.getType().toString();
    }


    /**                         1503292
     * Return the allele frequency as annotated with an AF or GMAF attribute.  A value of -1 indicates
     * no annotation (unknown allele frequency).
     */
    public double getAlleleFreq() {
        double alleleFreq = Double.parseDouble(variantContext.getAttributeAsString("AF", "-1"));
        if (alleleFreq < 0) {
            alleleFreq = Double.parseDouble(variantContext.getAttributeAsString("GMAF", "-1"));
        }
        return alleleFreq;

    }

    /**
     * Return the allele fraction for this variant.  The allele fraction is similiar to allele frequency, but is based
     * on the samples in this VCF as opposed to an AF or GMAF annotation.
     * <p/>
     * A value of -1 indicates unknown
     */
    public double getAlleleFraction() {

        int total = getHomVarCount() + getHetCount() + getHomRefCount();
        return total == 0 ? -1 : (((double) getHomVarCount() + ((double) getHetCount()) / 2) / total);
    }


    public Collection<String> getSampleNames() {
        return variantContext.getSampleNames();
    }

    public Map<String, Object> getAttributes() {
        return variantContext.getAttributes();
    }

    public Genotype getGenotype(String sample) {
        // TODO -- cache these rather than make a new object each call?
        return new VCFGenotype(variantContext.getGenotype(sample));
    }

    public Collection<String> getFilters() {
        return variantContext.getFilters();
    }

    public int getHomVarCount() {
        return zygosityCount.getHomVar();
    }

    public int getHetCount() {
        return zygosityCount.getHet();
    }

    public int getHomRefCount() {
        return zygosityCount.getHomRef();
    }

    public int getNoCallCount() {
        return zygosityCount.getNoCall();
    }

    public String getChr() {
        return chr;
    }

    public int getStart() {
        return isIndel ? variantContext.getStart(): variantContext.getStart() - 1;
    }

    public int getEnd() {
        return variantContext.getEnd();
    }

    @Override
    public String toString() {
        return String.format("VCFVariant[%s:%d-%d]",getChr(),getStart(),getEnd());
    }

    /**
* @author Jim Robinson
* @date Aug 1, 2011
*/
public static class ZygosityCount {
    private int homVar = 0;
    private int het = 0;
    private int homRef = 0;
    private int noCall = 0;

    public void incrementCount(Genotype genotype) {
        if (genotype != null) {
            if (genotype.isHomVar()) {
                homVar++;
            } else if (genotype.isHet()) {
                het++;
            } else if (genotype.isHomRef()) {
                homRef++;
            } else {
                noCall++;
            }
        }
    }

    public int getHomVar() {
        return homVar;
    }

    public int getHet() {
        return het;
    }

    public int getHomRef() {
        return homRef;
    }

    public int getNoCall() {
        return noCall;
    }

}
}
