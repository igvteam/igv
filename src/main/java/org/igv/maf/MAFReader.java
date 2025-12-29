/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.igv.maf;

import java.io.IOException;
import java.util.Collection;
import java.util.List;

/**
 * @author jrobinso
 */
public interface MAFReader {

    List<MultipleAlignmentBlock> loadAlignments(String chr, int start, int end) throws IOException;

    /**
     * Return the chromosome names represented in this file.   Can return null if unknown.
     *
     * @return
     */
    Collection<String> getChrNames();

    Collection<String> getSpecies();

    String getSpeciesName(String speciesId);

    String getRefId();
}
