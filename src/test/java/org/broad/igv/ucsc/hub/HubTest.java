package org.broad.igv.ucsc.hub;

import org.broad.igv.feature.genome.load.GenomeConfig;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;
import java.util.List;

import static org.junit.Assert.*;

public class HubTest {

    @Test
    public void testGetGenomeConfig() throws IOException {

        String hubFile = TestUtils.DATA_DIR + "hubs/hub.txt";
        Hub hub =  Hub.loadHub(hubFile);
        assertNotNull(hub.hubStanza);
        assertNotNull(hub.genomeStanza);
        assertEquals(22, hub.trackStanzas.size());

        GenomeConfig genomeConfig = hub.getGenomeConfig(true);
        assertNotNull(genomeConfig);
        assertEquals("GCF_000186305.1", genomeConfig.getId());
        assertEquals("Python bivittatus (GCF_000186305.1)", genomeConfig.getName());
        assertNotNull(genomeConfig.getTwoBitBptURL());
        assertNotNull(genomeConfig.getTwoBitURL());
        assertNotNull(genomeConfig.getChromAliasBbURL());
        assertNotNull(genomeConfig.getCytobandBbURL());
    }

    @Test
    public void testGetGroupedTrackConfigurations() throws IOException {
        String hubFile = TestUtils.DATA_DIR + "hubs/hub.txt";
        Hub hub =  Hub.loadHub(hubFile);
        List<TrackConfigGroup> groupedTrackConfigurations  = hub.getGroupedTrackConfigurations();
        assertEquals(5, groupedTrackConfigurations.size());
    }


    /**
     * track hgUnique
     * type bigBed 3
     * html hgUnique.html
     * visibility pack
     * shortLabel CHM13 unique
     * group map
     * compositeTrack on
     * priority 3
     * longLabel CHM13 unique in comparison to GRCh38/hg38 and GRCh37/hg19
     *
     * track hgUniqueHg38
     * type bigBed 3
     * shortLabel CHM13 unique for hg38
     * group compGeno
     * priority 1
     * bigDataUrl /gbdb/hs1/hgUnique/hgUnique.hg38.bb
     * parent hgUnique on
     * color 100,0,200
     * longLabel CHM13 unique in comparison to GRCh38/hg38
     *
     * track hgUniquehg19
     * type bigBed 3
     * shortLabel CHM13 unique for hg19
     * group compGeno
     * priority 2
     * bigDataUrl /gbdb/hs1/hgUnique/hgUnique.hg19.bb
     * parent hgUnique off
     * color 200,0,100
     * longLabel CHM13 unique in comparison to GRCh37/hg19
     *
     * @throws Exception
     */
    @Test
    public void testHS1() throws Exception {

        String hubFile = "https://hgdownload.soe.ucsc.edu/gbdb/hs1/hubs/public/hub.txt";

        Hub hub = Hub.loadHub(hubFile);
        List<TrackConfigGroup> groups = hub. getGroupedTrackConfigurations();

        // Find the compGeno group
        TrackConfigGroup mapGroup = null;
        for(TrackConfigGroup group : groups) {
            if(group.name.equals("map")) {
                mapGroup = group;
                break;
            }
        }
        assertNotNull(mapGroup);
        assertEquals(5, mapGroup.trackContainers.size());



    }
 }