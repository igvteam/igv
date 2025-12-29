package org.igv.feature.genome.load;

import org.junit.Test;

import static org.junit.Assert.*;

public class GenomeHubLoaderTest {

    @Test
    public void convert() {
        String gcfAccession = "GCF_000442705.1";
        String gcfURL = "https://hgdownload.soe.ucsc.edu/hubs/GCF/000/442/705/GCF_000442705.1/hub.txt";
        String hubTxt = HubGenomeLoader.convertToHubURL(gcfAccession);
        assertEquals(gcfURL, hubTxt);
    }

}