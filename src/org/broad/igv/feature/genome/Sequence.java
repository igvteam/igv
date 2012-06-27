package org.broad.igv.feature.genome;

import java.util.Collection;
import java.util.List;
import java.util.Set;

/**
 * @author jrobinso
 * @date 8/8/11
 *
 * Represents a reference sequence.
 */
public interface Sequence {

    byte[] getSequence(String chr, int start, int end);

    public byte getBase(String chr, int position);

    List<String> getChromosomeNames();

    int getChromosomeLength(String chrname);
}
                                                                    