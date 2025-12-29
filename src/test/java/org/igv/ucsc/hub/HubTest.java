package org.igv.ucsc.hub;

import org.igv.feature.genome.load.GenomeConfig;
import org.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;
import java.util.List;

import static org.junit.Assert.*;

public class HubTest {

    @Test
    public void testGetGenomeConfigs() throws IOException {

        String hubFile = TestUtils.DATA_DIR + "hubs/hub.txt";
        Hub hub = HubParser.loadAssemblyHub(hubFile);
        assertNotNull(hub.hubStanza);
        //assertNotNull(hub.genomeStanza);
       // assertEquals(22, hub.trackStanzas.size());

        GenomeConfig genomeConfig = hub.getGenomeConfigs().get(0);
        assertNotNull(genomeConfig);
        assertEquals("GCF_000186305.1", genomeConfig.id);
        assertEquals("Python (GCF_000186305.1)", genomeConfig.getName());
        assertNotNull(genomeConfig.twoBitBptURL);
        assertNotNull(genomeConfig.twoBitURL);
        assertNotNull(genomeConfig.chromAliasBbURL);
        assertNotNull(genomeConfig.cytobandBbURL);
    }

    @Test
    public void testGetGroupedTrackConfigurations() throws IOException {
        String hubFile = TestUtils.DATA_DIR + "hubs/hub.txt";
        Hub hub = HubParser.loadAssemblyHub(hubFile);
        List<TrackConfigContainer> groupedTrackConfigurations = hub.getGroupedTrackConfigurations("GCF_000186305.1");
        assertEquals(5, groupedTrackConfigurations.size());
    }

    @Test
    public void testNCBIHostedHub() throws IOException {

        String hubFile = "https://ftp.ncbi.nlm.nih.gov/snp/population_frequency/TrackHub/latest/hub.txt";
        Hub hub = HubParser.loadHub(hubFile);
        List<TrackConfigContainer> groupedTrackConfigurations = hub.getGroupedTrackConfigurations("hg38");
        assertEquals(1, groupedTrackConfigurations.size());
        assertEquals(12, groupedTrackConfigurations.get(0).tracks.size());

    }
}
