package org.igv.variant;

/**
 * Represents an allele.
 *
 * Note:  Currently interface could easily be replaced by a String, but we expect more complex alleles when the
 * structural variant code is added.
 *
 * @author Jim Robinson
 * @date Aug 1, 2011
 */
public interface Allele {
    public byte[] getBases() ;
    public String getDisplayString();
}
