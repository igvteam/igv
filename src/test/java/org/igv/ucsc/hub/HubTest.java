package org.igv.ucsc.hub;

import org.igv.feature.Chromosome;
import org.igv.feature.genome.Genome;
import org.igv.feature.genome.load.GenomeConfig;
import org.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.*;

public class HubTest {

    @Test
    public void testGetGenomeConfigs() throws IOException {

        String hubFile = TestUtils.DATA_DIR + "hubs/hub.txt";
        Hub hub = HubParser.loadHub(hubFile);
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
        Hub hub = HubParser.loadHub(hubFile);
        List<TrackConfigContainer> groupedTrackConfigurations = hub.getGroupedTrackConfigurations("GCF_000186305.1");
        assertEquals(5, groupedTrackConfigurations.size());
    }

    /**
     * A hub can mix UCSC hosted genomes referenced by id with custom assemblies defined by a twoBitPath.  Only the
     * latter yield genome configs.
     */
    @Test
    public void testMixedHub() throws IOException {
        String hubFile = TestUtils.DATA_DIR + "hubs/mixed/hub.txt";
        Hub hub = HubParser.loadHub(hubFile);

        assertEquals(Arrays.asList("hg38", "myAsm1"), hub.getGenomeIds());
        assertFalse(hub.hasTwoBit("hg38"));
        assertTrue(hub.hasTwoBit("myAsm1"));
        assertNull(hub.getGenomeConfig("hg38"));
        assertEquals(1, hub.getGenomeConfigs().size());

        GenomeConfig config = hub.getGenomeConfig("myAsm1");
        assertEquals("myAsm1", config.id);
        assertTrue(config.twoBitURL.endsWith("myAsm1/myAsm1.2bit"));
        assertEquals("Testus (myAsm1)", config.getName());   // HubParser keeps the first word of most values

        // An IGV genome is matched to the hub by its UCSC id, so GRC names and the 1kg genomes match an hg38 hub.
        Genome grch38 = new Genome("GRCh38", List.of(new Chromosome(0, "chr1", 1000)));
        assertEquals("hg38", grch38.getUCSCId());
        assertTrue(hub.referencesGenome(grch38.getUCSCId()));
        assertFalse(hub.referencesGenome("hg19"));
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
