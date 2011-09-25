package org.broad.igv.hic.tools;

import java.util.Iterator;

/**
 * @author Jim Robinson
 * @date 9/24/11
 */
public interface PairIterator extends Iterator<AlignmentPair> {
    boolean hasNext();

    AlignmentPair next();

    void remove();

    void close();
}
