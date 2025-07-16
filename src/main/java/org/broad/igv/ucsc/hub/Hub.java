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
    private List<GenomeConfig> genomeConfigs;

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

    public boolean isOneFile() {
        return "on".equals(hubStanza.getProperty("useOneFile"));
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

        TrackDbHub trackHub = getTrackDbHub(genomeId);
        if (trackHub == null) {
            log.error("Could not get track hub for genome " + genomeId);
            return Collections.emptyList();
        }
        String longLabel = this.getLongLabel();
        String hubLabel = longLabel != null && longLabel.length() < 50 ? longLabel : this.getShortLabel();
        return trackHub.getGroupedTrackConfigurations(hubLabel);
    }

    private TrackDbHub getTrackDbHub(String genomeId) {
        TrackDbHub trackHub = trackHubMap.get(genomeId);
        if (trackHub == null && idMappings.containsKey(genomeId)) {
            trackHub = trackHubMap.get(idMappings.get(genomeId));
        }

        if (trackHub == null) {
            String alias = idMappings.get(genomeId);
            for (Stanza s : genomeStanzas) {
                if (genomeId.equals(s.getProperty("genome")) || (alias != null && alias.equals(s.getProperty("genome")))) {
                    try {
                        String trackDbURL = s.getProperty("trackDb");
                        if (trackDbURL == null) {
                            log.error("No trackDb property found for genome " + genomeId);
                            return null;
                        }
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
     * Return the genome configuration for this assembly hub.
     *
     * @return GenomeConfig
     */
    public List<GenomeConfig> getGenomeConfigs() {

        if (genomeConfigs == null) {

            genomeConfigs = new ArrayList<>();

            for (Stanza genomeStanza : genomeStanzas) {

                GenomeConfig genomeConfig = new GenomeConfig();
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
                if (isOneFile()) {
                    TrackDbHub trackDbHub = getTrackDbHub(genomeConfig.id);
                    if (trackDbHub != null) {
                        genomeConfig.cytobandBbURL = trackDbHub.findCytobandURL();
                    }
                }

                genomeConfig.getName(); // This will set the name if it is not already set.  Bad use of side effect.

                genomeConfigs.add(genomeConfig);
            }
        }

        return genomeConfigs;
    }

    private static final Map<String, String> idMappings = Map.ofEntries(
            Map.entry("hg38", "GCF_000001405.40"),
            Map.entry("mm39", "GCF_000001635.27"),
            Map.entry("mm10", "GCF_000001635.26"),
            Map.entry("bosTau9", "GCF_002263795.1"),
            Map.entry("canFam4", "GCF_011100685.1"),
            Map.entry("canFam6", "GCF_000002285.5"),
            Map.entry("ce11", "GCF_000002985.6"),
            Map.entry("dm6", "GCF_000001215.4"),
            Map.entry("galGal6", "GCF_000002315.6"),
            Map.entry("gorGor6", "GCF_008122165.1"),
            Map.entry("macFas5", "GCA_000364345.1"),
            Map.entry("panTro6", "GCA_002880755.3"),
            Map.entry("rn6", "GCF_000001895.5"),
            Map.entry("rn7", "GCF_015227675.2"),
            Map.entry("sacCer3", "GCF_000146045.2"),
            Map.entry("sacCer2", "GCF_000146045.2"),
            Map.entry("susScr11", "GCF_000003025.6"),
            Map.entry("taeGut1", "GCF_000002275.3"),
            Map.entry("tetNig2", "GCF_000002275.3"),
            Map.entry("xenTro10", "GCF_000002035.6"),
            Map.entry("xenTro9", "GCF_000002035.6"),
            Map.entry("tair10", "GCF_000001735.4")
    );
}
