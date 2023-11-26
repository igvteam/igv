package org.broad.igv.feature.genome.load;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ucsc.Hub;

import java.io.IOException;

public class HubGenomeLoader extends GenomeLoader {

    String hubURL;

    public HubGenomeLoader(String hubURL) {
        this.hubURL = hubURL;
    }

    public static boolean isHubURL(String obj) {
        return obj.endsWith("/hub.txt");
    }

    @Override
    public Genome loadGenome() throws IOException {

        Hub hub = Hub.loadHub(this.hubURL);

        GenomeConfig config = hub.getGenomeConfig(null);

        return new Genome(config);

    }
}
