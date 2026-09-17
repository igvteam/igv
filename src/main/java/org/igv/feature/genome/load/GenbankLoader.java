package org.igv.feature.genome.load;

import org.igv.Globals;
import org.igv.feature.genome.Genome;
import org.igv.feature.genome.InMemorySequence;
import org.igv.feature.genome.Sequence;
import org.igv.track.FeatureTrack;
import org.igv.ui.IGV;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GenbankLoader extends GenomeLoader {

    private String genomePath;

    public GenbankLoader(String genomePath) {
        this.genomePath = genomePath;
    }

    @Override
    public Genome loadGenome() throws IOException {
       // Genome newGenome;

        GenbankParser genbankParser = new GenbankParser(genomePath);
        genbankParser.readFeatures(true);

        String name = genbankParser.getLocusName();
        String chr = genbankParser.getChr();
        if (!name.equals(chr)) {
            name = name + " (" + chr + ")";
        }
        byte[] seq = genbankParser.getSequence();
        Sequence sequence = new InMemorySequence(chr, seq);

        GenomeConfig config = new GenomeConfig();

        // The genome ID is the versioned accession ("NC_012920.1"), which is how the genome is listed and so how a
        // session or batch command will name it.  The sequence keeps the unversioned name.
        config.id = genbankParser.getGenomeId();
        config.setName(name);
        config.setSequence(sequence);

        String[] aliases = genbankParser.getAliases();
        if (aliases != null) {
            List<String> aliasList = new ArrayList<String>();
            aliasList.add(chr);
            for (String a : aliases) {
                aliasList.add(a);
            }
            config.setChromAliases((Arrays.asList(aliasList)));
        }


        Genome genome =  new Genome(config);

        if (IGV.hasInstance() && !Globals.isHeadless()) {
            genome.getFeatureDB().addFeatures(genbankParser.getFeatures());
            FeatureTrack geneFeatureTrack = createGeneTrack(genome, genbankParser.getFeatures());
            genome.setGeneTrack(geneFeatureTrack);
        }


        return genome;
    }
}
