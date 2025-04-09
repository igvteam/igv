package org.broad.igv.feature.genome.load;

import org.broad.igv.Globals;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.InMemorySequence;
import org.broad.igv.feature.genome.Sequence;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.ui.IGV;

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
        config.id = (chr);
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
            FeatureDB.addFeatures(genbankParser.getFeatures(), genome);
            FeatureTrack geneFeatureTrack = createGeneTrack(genome, genbankParser.getFeatures());
            genome.setGeneTrack(geneFeatureTrack);
        }


        return genome;
    }
}
