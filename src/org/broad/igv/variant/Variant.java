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

import org.broad.tribble.Feature;

import java.util.Collection;
import java.util.Map;
import java.util.Set;

/**
 * Represents a variant call on a collection of samples.  
 *
 * @author jrobinso
 * @date Jan 20, 2011
 */
public interface Variant extends Feature {

    /**
     * @return the identifier for this variant.  This is the "ID" attribute in VCF files, and is typically a
     * dbSnp id.   This method can return null, which means that no ID has been assigned.  There is not
     */
    String getID();

    /**
     * @return   the type of this variant as a String.  Any String value is legal, typical values from VCF files include
     *      NO_VARIATION, SNP, MNP, INDEL, SYMBOLIC, MIXED
     */
    String getType();


    /**
     *
     * @return return true if the variant has been marked as filtered
     */
    boolean isFiltered();

    /**
=
     * @return the phred scale quality score of this variant.
     */
    double getPhredScaledQual();


    /**

     * @return the map of all attributes.  This cannot be null, if there are no attributes return an empty map.
     */
    Map<String, Object> getAttributes();

    /**
     * Return the string value of an attribute.  Can be null.
     *
     * @param key
     * @return the attribute for the given key
     */
    String getAttributeAsString(String key);

    /**
     * @return The reference sequence for this variant
     */
    String getReference();

    /**
     * Return the set of alternate alleles for this variant.  This should not return null.
     *
     * @return The set of alternate alleles
     */
    Set<Allele> getAlternateAlleles();

    /**
     * Return the allele frequency for this variant, possibly from an annotation as opposed to the actual
     * sample genotypes associated with this variant.
     */
    double getAlleleFreq();

    /**
     * Return the allele fraction for this variant.  The allele fraction is similiar to allele frequency,
     * but is based on the specific samples in this VCF as opposed to an annotation.
     * <p/>
     * A value of -1 indicates unknown
     */
    double getAlleleFraction();

    /**
     * Return the list of sample names associated with this variant.
     *
     * @return
     */
    Collection<String> getSampleNames();

    /**
     * Return the genotype for the given sample
     *
     * @param sample
     * @return
     */
    Genotype getGenotype(String sample);

    /**
     * @return the list of filters applied to this variant.  Should not return null, return an empty collection if no filters.
     */
    Collection<String> getFilters();

    /**
     * @return the count of genotypes for this variant called as homozygous variant
     */
    public int getHomVarCount();

    /**
     * @return the count of genotypes for this variant called as heterozygous variant
     */
    public int getHetCount();

    /**
     * @return the count of genotypes for this variant called as homozygous reference
     */
    public int getHomRefCount();

    /**
     * @return the count of genotypes for this variant that are no-calls
     */
    public int getNoCallCount();

}
