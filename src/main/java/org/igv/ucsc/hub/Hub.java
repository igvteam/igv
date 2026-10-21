package org.igv.ucsc.hub;

import org.igv.feature.genome.load.GenomeConfig;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;

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
                this.hubStanza.getProperty("shortLabel"),
                this.hubStanza.getProperty("longLabel"),
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
     * Return the ids of the genomes referenced by this hub, in genomes.txt order.  Each is a UCSC database name
     * (hg38), a GenArk accession (GCF_...), or the id of a custom assembly defined by the hub itself.
     */
    public List<String> getGenomeIds() {
        return genomeStanzas.stream()
                .map(stanza -> stanza.getProperty("genome"))
                .filter(Objects::nonNull)
                .collect(Collectors.toList());
    }

    /**
     * Return true if the genome stanza for genomeId defines its own sequence (a "twoBitPath"), that is, the hub
     * defines the assembly rather than referencing a UCSC hosted one.
     */
    public boolean hasTwoBit(String genomeId) {
        Stanza stanza = getGenomeStanza(genomeId);
        return stanza != null && stanza.hasProperty("twoBitPath");
    }

    /**
     * Return true if this hub references the genome with the given UCSC id, either directly or by the accession
     * IGV maps it to (e.g. hg38 and GCF_000001405.40).  Callers pass Genome.getUCSCId(), never the IGV id, so
     * aliases such as hg38_1kg or GRCh38 match an hg38 hub.
     */
    public boolean referencesGenome(String ucscId) {
        return getGenomeStanza(ucscId) != null;
    }

    private Stanza getGenomeStanza(String genomeId) {
        if (genomeId == null) {
            return null;
        }
        String alias = idMappings.get(genomeId);
        for (Stanza s : genomeStanzas) {
            String g = s.getProperty("genome");
            if (genomeId.equals(g) || (alias != null && alias.equals(g))) {
                return s;
            }
        }
        return null;
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
            Stanza s = getGenomeStanza(genomeId);
            if (s != null) {
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
     * Return genome configurations for the assemblies defined by this hub, that is, the genome stanzas with a
     * "twoBitPath".  Stanzas that reference a UCSC hosted genome by id only have no sequence and are skipped.
     *
     * @return GenomeConfig
     */
    public List<GenomeConfig> getGenomeConfigs() {

        if (genomeConfigs == null) {

            genomeConfigs = new ArrayList<>();

            for (Stanza genomeStanza : genomeStanzas) {

                if (!genomeStanza.hasProperty("twoBitPath")) {
                    continue;
                }

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

    /**
     * Return the genome configuration for the assembly genomeId defined by this hub, or null if the hub does not
     * define that assembly (unknown id, or a UCSC hosted genome referenced by id only).
     */
    public GenomeConfig getGenomeConfig(String genomeId) {
        Stanza stanza = getGenomeStanza(genomeId);
        if (stanza == null) {
            return null;
        }
        String id = stanza.getProperty("genome");
        return getGenomeConfigs().stream().filter(c -> id.equals(c.id)).findFirst().orElse(null);
    }

    /**
     * UCSC database name -> GenArk accession, for the UCSC genomes IGV hosts that also have a GenArk hub.  Lets an
     * IGV genome such as "mm39" find a hub stanza named by accession.  Verified against UCSC GenArk Oct 2026;
     * hosted genomes without a GenArk hub (hg19, hg18, mm9, dm3, ...) are omitted.
     */
    private static final Map<String, String> idMappings = Map.ofEntries(
            Map.entry("hg38", "GCF_000001405.40"),
            Map.entry("hs1", "GCA_009914755.4"),
            Map.entry("mm39", "GCF_000001635.27"),
            Map.entry("mm10", "GCF_000001635.26"),
            Map.entry("rn6", "GCF_000001895.5"),
            Map.entry("rn7", "GCF_015227675.2"),
            Map.entry("danRer11", "GCF_000002035.6"),
            Map.entry("ce11", "GCF_000002985.6"),
            Map.entry("dm6", "GCF_000001215.4"),
            Map.entry("sacCer3", "GCF_000146045.2"),
            Map.entry("bosTau8", "GCF_000003055.6"),
            Map.entry("bosTau9", "GCF_002263795.1"),
            Map.entry("canFam3", "GCF_000002285.3"),
            Map.entry("canFam4", "GCF_011100685.1"),
            Map.entry("canFam6", "GCF_000002285.5"),
            Map.entry("panTro6", "GCA_002880755.3"),
            Map.entry("gorGor6", "GCF_008122165.1"),
            Map.entry("panPan2", "GCF_000258655.2"),
            Map.entry("macFas5", "GCA_000364345.1"),
            Map.entry("galGal6", "GCF_000002315.6"),
            Map.entry("susScr11", "GCF_000003025.6"),
            Map.entry("tair10", "GCF_000001735.4"),
            Map.entry("xenTro10", "GCF_000004195.4"),
            Map.entry("tetNig2", "GCA_000180735.1")
    );
}
