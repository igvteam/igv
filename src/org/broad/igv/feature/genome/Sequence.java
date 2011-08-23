package org.broad.igv.feature.genome;

/**
 * @author jrobinso
 * @date 8/8/11
 *
 * Represents a reference sequence.
 */
public interface Sequence {

    byte[] readSequence(String chr, int start, int end);
}
