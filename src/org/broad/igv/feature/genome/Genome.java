/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.feature.genome;

import org.broad.igv.dev.api.api;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.track.FeatureTrack;

import java.util.Collection;
import java.util.List;

/**
 * @author Jim Robinson
 * @date 1/21/12
 */
public interface Genome {

    String getId();

    String getSpecies();

    String getHomeChromosome();

    Chromosome getChromosome(String chrName);

    Collection<Chromosome> getChromosomes();

    List<String> getAllChromosomeNames();

    String getNextChrName(String chr);

    String getPrevChrName(String chr);

    String getChromosomeAlias(String str);

    long getTotalLength();

    /**
     * Return the nucleotide sequence on the + strand for the genomic interval.  This method can return null
     * if sequence is not available.
     *
     * @param chr
     * @param start  start position in "zero-based" coordinates
     * @param end  end position
     * @return  sequence, or null if not available
     */
    @api
    byte[] getSequence(String chr, int start, int end);

    String getDisplayName();

    byte getReference(String chr, int pos);

    FeatureTrack getGeneTrack();


    // Methods to support whole-genome view follow

    /**
     * Return "getChromosomeNames()" with small chromosomes removed.
     *
     * @return
     */
    List<String> getLongChromosomeNames();

    long getNominalLength();

    long getCumulativeOffset(String chr);

    /**
     * Covert the chromosome coordinate in basepairs to genome coordinates in kilo-basepairs
     *
     * @param chr
     * @param locationBP
     * @return The overall genome coordinate, in kilo-bp
     */
    int getGenomeCoordinate(String chr, int locationBP);

    /**
     * Translated a genome coordinate, in kilo-basepairs, to a chromosome & position in basepairs.
     *
     * @param genomeKBP The "genome coordinate" in kilo-basepairs.  This is the distance in kbp from the start of the
     *                  first chromosome.
     * @return
     */
    ChromosomeCoordinate getChromosomeCoordinate(int genomeKBP);
}
