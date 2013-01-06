package org.broad.igv.feature.genome;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.Cytoband;
import org.broad.igv.track.FeatureTrack;

import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Jim Robinson
 * @date 1/21/12
 */
public interface Genome {

    String getId();

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
    long getNominalLength();

    List<String> getLongChromosomeNames();

    long getCumulativeOffset(String chr);

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
