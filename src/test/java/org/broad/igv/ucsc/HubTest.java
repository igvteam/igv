package org.broad.igv.ucsc;

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
        assertNotNull(hub.hub);
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

}