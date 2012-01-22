package org.broad.igv.feature.genome;

import org.broad.igv.feature.Chromosome;

import java.util.Collection;
import java.util.List;

/**
 * @author Jim Robinson
 * @date 1/21/12
 */
public interface Genome {

    String getId();

    String getHomeChromosome();

    Chromosome getChromosome(String chrName);

    Collection<Chromosome> getChromosomes();

    List<String> getChromosomeNames();

    String getChromosomeAlias(String str);

    long getLength();

    long getCumulativeOffset(String chr);

    int getGenomeCoordinate(String chr, int locationBP);

    /**
     * Translated a genome coordinate, in kilo-basepairs, to a chromosome & position in basepairs
     *
     * @param genomeKBP The "genome coordinate" in kilo-basepairs.  This is the distance in kbp from the start of the
     *                  first chromosome.
     * @return
     */
    ChromosomeCoordinate getChromosomeCoordinate(int genomeKBP);

    String getNextChrName(String chr);

    String getPrevChrName(String chr);

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
}
