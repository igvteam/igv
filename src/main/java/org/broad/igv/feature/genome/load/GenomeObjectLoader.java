package org.broad.igv.feature.genome.load;

import org.broad.igv.feature.CytoBandFileParser;
import org.broad.igv.feature.genome.*;
import org.broad.igv.feature.genome.fasta.FastaBlockCompressedSequence;
import org.broad.igv.feature.genome.fasta.FastaIndexedSequence;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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
