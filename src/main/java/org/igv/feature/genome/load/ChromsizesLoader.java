package org.igv.feature.genome.load;

import org.igv.feature.Chromosome;
import org.igv.feature.genome.Genome;

import java.io.IOException;
import java.util.List;

public class ChromsizesLoader extends GenomeLoader {

    private String genomePath;

    public ChromsizesLoader(String genomePath) {
        this.genomePath = genomePath;
    }

    /**
     * Define a minimal genome from a chrom.sizes file.
     *
     * @return
     * @throws IOException
     */
    @Override
    public Genome loadGenome() throws IOException {
        String genomeId = genomePath;
        List<Chromosome> chromosomes = ChromSizesParser.parse(genomePath);
        Genome newGenome = new Genome(genomeId, chromosomes);
        return newGenome;

    }
}
