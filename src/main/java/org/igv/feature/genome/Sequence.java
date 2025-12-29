package org.igv.feature.genome;

import org.igv.feature.Chromosome;
import org.igv.ui.panel.ReferenceFrame;

import java.util.List;

/**
 * @author jrobinso
 * @date 8/8/11
 *
 * Represents a reference sequence.
 */
public interface Sequence {

    /**
     * Return the sequence for the given range.  If sequence named "seq" does not exist returns null.
     * @param seq
     * @param start
     * @param end
     * @return  The sequence in bytes, or null if no sequence exists
     */
    byte[] getSequence(String seq, int start, int end) ;

    byte getBase(String seq, int position);

    List<String> getChromosomeNames();

    /**
     * Return the given sequence length.  If no sequence exists with name "seq" return -1.
     *
     * @param seq
     * @return
     */
    int getChromosomeLength(String seq);

    List<Chromosome> getChromosomes();

    default boolean hasChromosomes() {
        return true;
    }

}
