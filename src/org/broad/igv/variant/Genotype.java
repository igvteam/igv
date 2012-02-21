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

package org.broad.igv.variant;

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
     *   NO_CALL, HOM_REF, HET, HOM_VAR, UNAVAILABLE.   Can be null.
     */
    String getType();

    /**
     * @return a string representation of this genotype.  Can be null.
     */
    String getGenotypeString();

    /**
     * @return a key-value map of this genotypes attributes.  Cannot be null.
     */
    Map<String, Object> getAttributes();

    /**
     * Return the attribute as a string.  Can be null.
     *
     * @param key
     * @return attribute
     */
    String getAttributeAsString(String key);

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
     * @return  the Phred scale quality score for this genotype
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
