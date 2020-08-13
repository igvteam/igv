package org.broad.igv.feature.genome.load;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.Genome;

import java.io.IOException;
import java.util.List;

public class ChromsizesLoader extends GenomeLoader {

    private String genomePath;

    public ChromsizesLoader(String genomePath) {
        this.genomePath = genomePath;
    }

    /**
     * Define a minimal genome from a chrom.sizes file.  It is assumed (required) that the file follow the
     * UCSC naming convention  =>  [id].chrom.sizes
     * @return
     * @throws IOException
     */
    @Override
    public Genome loadGenome() throws IOException {
        int firstPeriodIdx = genomePath.indexOf('.');
        String genomeId = genomePath.substring(0, firstPeriodIdx);
        List<Chromosome> chromosomes = ChromSizesParser.parse(genomePath);
        Genome newGenome = new Genome(genomeId, chromosomes);
        return newGenome;

    }
}
