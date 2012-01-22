package org.broad.igv.feature.genome;

import org.broad.igv.feature.Chromosome;

import java.util.Collection;
import java.util.List;

/**
 * @author Jim Robinson
 * @date 1/21/12
 */
public interface Genome {
    String getChromosomeAlias(String str);

    String getHomeChromosome();

    Chromosome getChromosome(String chrName);

    List<String> getChromosomeNames();

    Collection<Chromosome> getChromosomes();

    long getLength();

    long getCumulativeOffset(String chr);

    int getGenomeCoordinate(String chr, int locationBP);

    GenomeImpl.ChromosomeCoordinate getChromosomeCoordinate(int genomeKBP);

    String getId();

    String getNextChrName(String chr);

    String getPrevChrName(String chr);

    byte[] getSequence(String chr, int start, int end);

    String getDisplayName();

    byte getReference(String chr, int pos);
}
