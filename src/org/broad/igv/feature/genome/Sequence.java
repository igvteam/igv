package org.broad.igv.feature.genome;

/**
 * @author jrobinso
 * @date 8/8/11
 *
 * Represents a reference sequence.
 */
public interface Sequence {

    byte[] getSequence(String chr, int start, int end);

    public byte getBase(String chr, int position);

}
                                                                    