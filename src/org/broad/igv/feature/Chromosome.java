package org.broad.igv.feature;

import java.util.List;

/**
 * @author Jim Robinson
 * @date 1/21/12
 */
public interface Chromosome {
    int getLength();

    String getName();

    List<Cytoband> getCytobands();
}
