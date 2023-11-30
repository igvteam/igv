package org.broad.igv.feature.genome.load;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;

import java.io.IOException;

public class GenomeObjectLoader extends GenomeLoader {

    private static Logger log = LogManager.getLogger(GenomeObjectLoader.class);

    private GenomeConfig genomeConfig;

    public GenomeObjectLoader(GenomeConfig genomeConfig) throws IOException {
        this.genomeConfig = genomeConfig;
    }

    @Override
    public Genome loadGenome() throws IOException {
        return new Genome(genomeConfig);
    }

}
