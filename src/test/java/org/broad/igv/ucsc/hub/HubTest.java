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
        Hub hub = HubParser.loadAssemblyHub(hubFile);
        assertNotNull(hub.hubStanza);
        assertNotNull(hub.genomeStanza);
        assertEquals(22, hub.trackStanzas.size());

        GenomeConfig genomeConfig = hub.getGenomeConfig();
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
        List<TrackConfigContainer> groupedTrackConfigurations = hub.getGroupedTrackConfigurations();
        assertEquals(5, groupedTrackConfigurations.size());
    }

    @Test
    public void testNCBIHostedHub() throws IOException {

        String hubFile = "https://ftp.ncbi.nlm.nih.gov/snp/population_frequency/TrackHub/latest/hub.txt";
        Hub hub = HubParser.loadHub(hubFile, "hg38");
        List<TrackConfigContainer> groupedTrackConfigurations = hub.getGroupedTrackConfigurations();
        assertEquals(1, groupedTrackConfigurations.size());
        assertEquals(12, groupedTrackConfigurations.get(0).tracks.size());

    }
}
