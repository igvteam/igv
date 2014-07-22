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
