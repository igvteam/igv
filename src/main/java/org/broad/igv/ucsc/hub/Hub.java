package org.broad.igv.ucsc.hub;

import org.broad.igv.feature.genome.load.GenomeConfig;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Hub object representing a UCSC hub.  This is a collection of genome and track
 *
 *
 * Example genomes stanza or genomes.txt file for an assembly hub:
 *
 * genome ricCom1
 * trackDb ricCom1/trackDb.txt
 * groups ricCom1/groups.txt
 * description July 2011 Castor bean
 * twoBitPath ricCom1/ricCom1.2bit
 * organism Ricinus communis
 * defaultPos EQ973772:1000000-2000000
 * orderKey 4800
 * scientificName Ricinus communis
 * htmlPath ricCom1/description.html
 * transBlat yourLab.yourInstitution.edu 17777
 * blat yourLab.yourInstitution.edu 17779
 * isPcr yourLab.yourInstitution.edu 17779
 */


public class Hub {

    private static Logger log = LogManager.getLogger(Hub.class);
    private final String url;
    private int order = 0;
    Stanza hubStanza;
    List<Stanza> genomeStanzas;
    List<Stanza> groupStanzas;
    Map<String, TrackDbHub> trackHubMap;
    private GenomeConfig genomeConfig;

    Hub(String url,
        Stanza hubStanza,
        List<Stanza> genomeStanzas,
        List<Stanza> trackStanzas,
        List<Stanza> groupStanzas) throws IOException {

        this.url = url;
        this.hubStanza = hubStanza;
        this.genomeStanzas = genomeStanzas;
        this.groupStanzas = groupStanzas;
        this.trackHubMap = new HashMap<>();

        // trackStanzas will not be null if this is a "onefile" hub
        if (trackStanzas != null) {
            String genomeId = genomeStanzas.get(0).getProperty("genome");     // Assmuption here this is a single genome hub
            trackHubMap.put(genomeId, new TrackDbHub(trackStanzas, groupStanzas));
        }
    }

    public HubDescriptor getDescriptor() {
        String dbList = genomeStanzas.stream()
                .map(stanza -> stanza.getProperty("genome"))
                .collect(Collectors.joining(","));

        return new HubDescriptor(
                this.url,
                this.hubStanza.getProperty("longLabel"),
                this.hubStanza.getProperty("shortLabel"),
                dbList,
                this.hubStanza.getProperty("descriptionUrl")
        );
    }

    public int getOrder() {
        return order;
    }

    public void setOrder(int order) {
        this.order = order;
    }

    /**
     *  * Example genomes stanza or genomes.txt file for an assembly hub:
     *  *
     *  * genome ricCom1
     *  * trackDb ricCom1/trackDb.txt
     *  * groups ricCom1/groups.txt
     *  * description July 2011 Castor bean
     *  * twoBitPath ricCom1/ricCom1.2bit
     *  * organism Ricinus communis
     *  * defaultPos EQ973772:1000000-2000000
     *  * orderKey 4800
     *  * scientificName Ricinus communis
     *  * htmlPath ricCom1/description.html
     *  * transBlat yourLab.yourInstitution.edu 17777
     *  * blat yourLab.yourInstitution.edu 17779
     *  * isPcr yourLab.yourInstitution.edu 17779
     */

    public boolean isAssemblyHub() {
        return genomeStanzas.stream().allMatch(gs -> gs.hasProperty("twoBitPath"));
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

    public String getUrl() {
        return url;
    }

    public String getName() {
        return this.hubStanza.getProperty("shortLabel");
    }

    public String getDescription() {
        return this.hubStanza.getProperty("longLabel");
    }

    public List<TrackConfigContainer> getGroupedTrackConfigurations(String genomeId) {

        final TrackDbHub trackHub = getTrackDbHub(genomeId);
        String longLabel = this.getLongLabel();
        String hubLabel = longLabel != null && longLabel.length() < 50 ? longLabel : this.getShortLabel();
        return trackHub.getGroupedTrackConfigurations(hubLabel);
    }

    private TrackDbHub getTrackDbHub(String genomeId) {
        TrackDbHub trackHub = trackHubMap.get(genomeId);
        if (trackHub == null) {
            for (Stanza s : genomeStanzas) {
                if (genomeId.equals(s.getProperty("genome"))) {
                    try {
                        String trackDbURL = s.getProperty("trackDb");
                        List<Stanza> trackStanzas = HubParser.loadStanzas(trackDbURL);
                        trackHub = new TrackDbHub(trackStanzas, this.groupStanzas);
                        trackHubMap.put(genomeId, trackHub);
                    } catch (IOException e) {
                        log.error("Error loading trackDb file: " + s.getProperty("trackDb"), e);
                    }
                    //TODO -- groups
                }

            }
        }
        return trackHub;
    }

    public int getSupportedTrackCount(String genomeId) {
        TrackDbHub trackHub = getTrackDbHub(genomeId);
        if (trackHub == null) {
            return 0;
        }
        return trackHub.getSupportedTrackCount();
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


    /**
     * Return the genome configuration for this assembly hub.  Currently we support only a single genome
     * in an assembly hub.
     *
     * @return GenomeConfig
     */
    public GenomeConfig getGenomeConfig() {

        if (genomeConfig == null) {

            Stanza genomeStanza = genomeStanzas.get(0);
            if (genomeStanza == null) {
                throw new RuntimeException("No genome stanza found in hub");
            }
            genomeConfig = new GenomeConfig();
            genomeConfig.id = genomeStanza.getProperty("genome");
            genomeConfig.nameSet = "ucsc";
            genomeConfig.wholeGenomeView = false;
            genomeConfig.accession = genomeStanza.getProperty("genome");
            genomeConfig.taxId = genomeStanza.getProperty("taxId");
            genomeConfig.scientificName = genomeStanza.getProperty("scientificName");
            genomeConfig.twoBitURL = (genomeStanza.getProperty("twoBitPath"));
            genomeConfig.defaultPos = (genomeStanza.getProperty("defaultPos"));
            genomeConfig.blat = genomeStanza.getProperty("blat");
            genomeConfig.chromAliasBbURL = genomeStanza.getProperty("chromAliasBb");
            genomeConfig.twoBitBptURL = genomeStanza.getProperty("twoBitBptURL");
            if (genomeConfig.twoBitBptURL == null) {
                genomeConfig.twoBitBptURL = genomeStanza.getProperty("twoBitBptUrl");
            }
            genomeConfig.description = genomeStanza.getProperty("description");
            genomeConfig.organism = genomeStanza.getProperty("organism");
            genomeConfig.scientificName = genomeStanza.getProperty("scientificName");
            genomeConfig.infoURL = (genomeStanza.getProperty("htmlPath"));
            genomeConfig.chromSizesURL = (genomeStanza.getProperty("chromSizes"));


            // Search for cytoband
            TrackDbHub trackDbHub = getTrackDbHub(genomeConfig.id);
            genomeConfig.cytobandBbURL = trackDbHub.findCytobandURL();
        }

        return genomeConfig;
    }
}
