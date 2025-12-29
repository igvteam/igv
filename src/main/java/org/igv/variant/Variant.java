package org.igv.variant;

import htsjdk.tribble.Feature;

import java.util.Collection;
import java.util.List;
import java.util.Map;

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

    boolean hasLog10PError();

    /**
     * @return the type of this variant as a String.  Any String value is legal, typical values from VCF files include
     * NO_VARIATION, SNP, MNP, INDEL, SYMBOLIC, MIXED
     */
    String getType();

    /**
     * @return true if all alleles are NON_REF.  This is common with gvcf files
     * @return
     */
    boolean isNonRef();


    /**
     * @return return true if the variant has been marked as filtered
     */
    boolean isFiltered();

    /**
     * =
     *
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
     * Return the list of alternate alleles for this variant.  This should not return null.
     *
     * @return
     */
    List<Allele> getAlternateAlleles();

    /**
     * Return the allele frequency for this variant, possibly from an annotation as opposed to the actual
     * sample genotypes associated with this variant.
     */
    double[] getAlleleFreqs();

    /**
     * Return the allele fraction for this variant.  The allele fraction is similiar to allele frequency,
     * but is based on the specific samples in this VCF as opposed to an annotation.
     * <p/>
     * A value of -1 indicates unknown
     */
    //  double getAlleleFraction();

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

//    /**
//     * @return the count of genotypes for this variant called as homozygous variant
//     */
//    public int getHomVarCount();
//
//    /**
//     * @return the count of genotypes for this variant called as heterozygous variant
//     */
//    public int getHetCount();
//
//    /**
//     * @return the count of genotypes for this variant called as homozygous reference
//     */
//    public int getHomRefCount();
//
//    /**
//     * @return the count of genotypes for this variant that are no-calls
//     */
//    public int getNoCallCount();

    double getMethlationRate();

    double getCoveredSampleFraction();

    String getPositionString();

    int[] getAlleleCounts();

    double getAlternateAlleleFrequency();

    int getTotalAlleleCount();

    double getAlleleFraction();
}
