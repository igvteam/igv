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

package org.broad.igv.variant;

import htsjdk.variant.variantcontext.GenotypeType;

import java.util.List;
import java.util.Map;

/**
 * Represents a genotype,  that is a variant call on a single sample
 *
 * @author Jim Robinson
 * @date Aug 1, 2011
 */
public interface Genotype {

    /**
     * @return the type.  Any string is legal,  typical values from VCF files are
     *         NO_CALL, HOM_REF, HET, HOM_VAR, UNAVAILABLE.   Can be null.
     */
    String getTypeString();

    /**
     * @return The GenotypeType.
     *         Should be GenotypeType.UNAVAILABLE if not known, never null
     *         See {@link htsjdk.variant.variantcontext.GenotypeType}
     */
    GenotypeType getType();

    /**
     * @return a string representation of this genotype.  Can be null.
     */
    String getGenotypeString();

    /**
     * @return a key-value map of this genotypes attributes.  Cannot be null.
     */
    Map<String, Object> getAttributes();


    /**
     * Return the attribute as a double,  or NaN if the attribute does not exist or cannot be converted to
     * a double.
     *
     * @param key
     * @return attribute
     */
    double getAttributeAsDouble(String key);

    /**
     * @return the alleles for this genotype.  Cannot be null, return an empty list if no alleles.
     */
    List<Allele> getAlleles();

    /**
     * @return the Phred scale quality score for this genotype
     */
    double getPhredScaledQual();

    /**
     * @return true if this genotype is homozygous variant
     */
    boolean isHomVar();

    /**
     * @return true if this genotype is heterozygous variant
     */
    boolean isHet();

    /**
     * @return true if this genotype is homozygous reference
     */
    boolean isHomRef();

    /**
     * @return true if this genotype is a no-call
     */
    boolean isNoCall();

}
