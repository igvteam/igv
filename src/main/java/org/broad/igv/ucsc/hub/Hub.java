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
    private int order = 0;
    Stanza hubStanza;
    Stanza genomeStanza;
    List<Stanza> trackStanzas;
    List<Stanza> groupStanzas;
    TrackDbHub trackHub;
    private GenomeConfig genomeConfig;

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

        if(genomeConfig == null) {

            genomeConfig = new GenomeConfig();
            genomeConfig.nameSet = "ucsc";
            genomeConfig.wholeGenomeView = false;
            genomeConfig.id = this.genomeStanza.getProperty("genome");
            genomeConfig.accession = this.genomeStanza.getProperty("genome");
            genomeConfig.taxId = this.genomeStanza.getProperty("taxId");
            genomeConfig.scientificName = this.genomeStanza.getProperty("scientificName");
            genomeConfig.twoBitURL = (this.genomeStanza.getProperty("twoBitPath"));
            genomeConfig.defaultPos = (this.genomeStanza.getProperty("defaultPos"));
            genomeConfig.blat = genomeStanza.getProperty("blat");
            genomeConfig.chromAliasBbURL = genomeStanza.getProperty("chromAliasBb");
            genomeConfig.twoBitBptURL = genomeStanza.getProperty("twoBitBptUrl");
            genomeConfig.description = genomeStanza.getProperty("description");
            genomeConfig.organism = genomeStanza.getProperty("organism");
            genomeConfig.scientificName = genomeStanza.getProperty("scientificName");
            genomeConfig.infoURL = (genomeStanza.getProperty("htmlPath"));
            genomeConfig.chromSizesURL = (genomeStanza.getProperty("chromSizes"));


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
                    genomeConfig.cytobandBbURL = t.getProperty("bigDataUrl");
                    break;
                }
            }
        }

        return genomeConfig;
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
