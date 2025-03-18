package org.broad.igv.ucsc.hub;

import org.broad.igv.feature.genome.load.GenomeConfig;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;

import java.io.IOException;
import java.util.*;


public class Hub {

    private static Logger log = LogManager.getLogger(Hub.class);
    private final String url;
    private final String trackDbURL;
    private int order;
    Stanza hubStanza;
    Stanza genomeStanza;
    List<Stanza> trackStanzas;
    List<Stanza> groupStanzas;
    TrackDbHub trackHub;

    Hub(String url,
        String trackDbURL,
        Stanza hubStanza,
        Stanza genomeStanza,
        List<Stanza> trackStanzas,
        List<Stanza> groupStanzas) throws IOException {

        this.url = url;
        this.trackDbURL = trackDbURL;

        // Validation checks
        this.hubStanza = hubStanza;
        this.genomeStanza = genomeStanza;
        this.trackStanzas = trackStanzas;
        this.groupStanzas = groupStanzas;

        // load includes.  Nested includes (includes within includes) are not supported
//        List<Stanza> includes = stanzas.stream().filter(s -> "include".equals(s.type)).toList();
//        for (Stanza s : includes) {
//            List<Stanza> includeStanzas = HubParser.loadStanzas(getDataURL(s.getProperty("include")));
//            for (Stanza inc : includeStanzas) {
//                inc.setProperty("visibility", "hide");
//                stanzas.add(inc);
//            }
//        }

        if (trackStanzas != null) {
            this.trackHub = new TrackDbHub(trackStanzas, groupStanzas);
        }
    }

    public int getOrder() {
        return order;
    }

    public void setOrder(int order) {
        this.order = order;
    }

    public boolean isAssemblyHub() {
        return genomeStanza.hasProperty("twoBitPath");
    }

    public String getShortLabel() {
        return hubStanza.hasProperty("shortLabel") ? hubStanza.getProperty("shortLabel") : this.url;
    }

    public String getLongLabel() {
        return hubStanza.hasProperty("longLabel") ? hubStanza.getProperty("longLabel") : this.url;
    }

    public String getDescriptionURL() {
        return hubStanza.hasProperty("descriptionUrl") ? hubStanza.getProperty("descriptionUrl") :
                hubStanza.hasProperty("desriptionUrl") ? hubStanza.getProperty("desriptionUrl") : null;
    }


    /**
     * Return the priority for the group.  The priority format is uncertain, but extends to at least 2 levels (e.g. 3.4).
     * Ignore levels > 3
     *
     * @param priorityString Priority as a string (e.g. 3.4)
     * @return A priority as an integer
     */

    static int getPriority(String priorityString) {
        try {
            String[] tokens = priorityString.trim().split("\\.");
            int p = Integer.parseInt(tokens[0]) * 100;
            if (tokens.length > 1) {
                p += Integer.parseInt(tokens[1]) * 10;
            }
            if (tokens.length > 2) {
                p += Integer.parseInt(tokens[2]);
            }
            return p;
        } catch (Exception e) {
            log.error("Error parsing priority string: " + priorityString, e);
            return Integer.MAX_VALUE;
        }
    }

    public GenomeConfig getGenomeConfig() {

        // TODO -- add blat?  htmlPath?

        GenomeConfig config = new GenomeConfig();
        config.setId(this.genomeStanza.getProperty("genome"));
        if (this.genomeStanza.hasProperty("scientificName")) {
            config.setName(this.genomeStanza.getProperty("scientificName"));
        } else if (this.genomeStanza.hasProperty("organism")) {
            config.setName(this.genomeStanza.getProperty("organism"));
        } else if (this.genomeStanza.hasProperty("description")) {
            config.setName(this.genomeStanza.getProperty("description"));
        }
        if (config.getName() == null) {
            config.setName(config.getId());
        } else {
            config.setName(config.getName() + " (" + config.getId() + ")");
        }

        config.setTwoBitURL(this.genomeStanza.getProperty("twoBitPath"));
        config.setNameSet("ucsc");
        config.setWholeGenomeView(false);

        if (this.genomeStanza.hasProperty("defaultPos")) {
            config.setDefaultPos(this.genomeStanza.getProperty("defaultPos"));
        }

        config.setDescription(config.getId());

        if (genomeStanza.hasProperty("blat")) {
            config.setBlat(genomeStanza.getProperty("blat"));
        }
        if (genomeStanza.hasProperty("chromAliasBb")) {
            config.setChromAliasBbURL(genomeStanza.getProperty("chromAliasBb"));
        }
        if (genomeStanza.hasProperty("twoBitBptURL")) {
            config.setTwoBitBptURL(genomeStanza.getProperty("twoBitBptURL"));
        }
        if (genomeStanza.hasProperty("twoBitBptUrl")) {
            config.setTwoBitBptURL(genomeStanza.getProperty("twoBitBptUrl"));
        }

        // chromSizes can take a very long time to load, and is not useful with the default WGV = off
        // if (this.genomeStanza.hasProperty("chromSizes")) {
        //     config.chromSizes = this.baseURL + this.genomeStanza.getProperty("chromSizes")
        // }

        if (genomeStanza.hasProperty("description")) {
            config.setDescription(config.getDescription() + "\n" + genomeStanza.getProperty("description"));
        }
        if (genomeStanza.hasProperty("organism")) {
            config.setDescription(config.getDescription() + "\n" + genomeStanza.getProperty("organism"));
        }
        if (genomeStanza.hasProperty("scientificName")) {
            config.setDescription(config.getDescription() + "\n" + genomeStanza.getProperty("scientificName"));
        }

        if (genomeStanza.hasProperty("htmlPath")) {
            config.setInfoURL(genomeStanza.getProperty("htmlPath"));
        }
        if (genomeStanza.hasProperty("chromSizes")) {
            config.setChromSizesURL(genomeStanza.getProperty("chromSizes"));
        }

        // Search for cytoband
        /*
        track cytoBandIdeo
        shortLabel Chromosome Band (Ideogram)
        longLabel Ideogram for Orientation
        group map
        visibility dense
        type bigBed 4 +
        bigDataUrl bbi/GCA_004027145.1_DauMad_v1_BIUU.cytoBand.bb
         */
        for (Stanza t : this.trackStanzas) {
            if ("cytoBandIdeo".equals(t.name) && t.hasProperty("bigDataUrl")) {
                config.setCytobandBbURL(t.getProperty("bigDataUrl"));
                break;
            }
        }

        return config;
    }

    public List<TrackConfigContainer> getGroupedTrackConfigurations() {
        if (this.trackHub == null) {
            try {
                List<Stanza> trackStanzas = HubParser.loadStanzas(this.trackDbURL);
                this.trackHub = new TrackDbHub(trackStanzas, this.groupStanzas);
            } catch (IOException e) {
                throw new RuntimeException("Error loading track configurations: " + e.getMessage(), e);
            }
        }
        String longLabel = this.getLongLabel();
        String hubLabel = longLabel != null && longLabel.length() < 50 ? longLabel : this.getShortLabel();
        return trackHub.getGroupedTrackConfigurations(hubLabel);
    }

    public String getUrl() {
        return url;
    }
}
