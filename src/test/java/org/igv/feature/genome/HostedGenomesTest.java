package org.igv.feature.genome;

import org.igv.feature.genome.load.GenomeConfig;
import org.junit.Test;
import static org.junit.Assert.*;

public class HostedGenomesTest {

    /**
     * A hosted genome is recognized by ID *and* sequence URL, so a user genome that reuses a hosted ID is not
     * mistaken for the hosted one.  Requires network access to the genome server.
     */
    @Test
    public void testIsIGVHosted() throws Exception {
        String hosted = org.igv.util.FileUtils.getContents(
            "https://raw.githubusercontent.com/igvteam/igv-data/refs/heads/main/genomes/json/hg38.json");
        assertTrue(HostedGenomes.isIGVHosted(GenomeConfig.fromJson(hosted)));

        // Same ID, different sequence -- a user genome reusing "hg38"
        String imposter = hosted.replace("https://igv.org/genomes/data/hg38/hg38.2bit",
                                         "https://example.org/mine/hg38.2bit");
        assertNotEquals(hosted, imposter);
        assertFalse(HostedGenomes.isIGVHosted(GenomeConfig.fromJson(imposter)));

        // Downloaded genome: twoBitURL rewritten to a local relative path
        String downloaded = hosted.replace("https://igv.org/genomes/data/hg38/hg38.2bit", "hg38/hg38.2bit");
        assertFalse(HostedGenomes.isIGVHosted(GenomeConfig.fromJson(downloaded)));

        // Unknown id
        String unknown = hosted.replace("\"id\": \"hg38\"", "\"id\": \"my_genome\"");
        assertFalse(HostedGenomes.isIGVHosted(GenomeConfig.fromJson(unknown)));

        assertFalse(HostedGenomes.isIGVHosted(null));
    }
}
