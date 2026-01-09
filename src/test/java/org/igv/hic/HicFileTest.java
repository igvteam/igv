package org.igv.hic;

import org.igv.feature.genome.Genome;
import org.igv.util.TestUtils;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.IOException;
import java.util.List;

import static org.junit.Assert.assertEquals;

public class HicFileTest {

    static String testURL = "https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/cc0accc9-a89f-4b40-bac6-de754a7f77c9/4DNFILI6XOB5.hic";
    static HicFile hicFile;

    @BeforeClass
    public static void setUp() throws Exception {
        Genome genome = TestUtils.mockUCSCGenome();
        hicFile = new HicFile(testURL, genome);
    }


    @Test
    public void getVersion() {
        int version = hicFile.getVersion();
        assertEquals(8, version);
    }

    @Test
    public void getContactRecords() throws IOException {
        Region region1 = new Region("chr1", 0, 249000000);
        Region region2 = new Region("chr1", 0, 249000000);
        int binSize = 250000;
        List<ContactRecord> records = hicFile.getContactRecords(region1, region2, "BP", binSize, "NONE", true, 0);
        assertEquals(396963, records.size());
    }

    @Test
    public void getNormalizationTypes() throws IOException {
        List<String> normTypes = hicFile.getNormalizationTypes();
        assertEquals(4, normTypes.size());
        assertEquals("NONE", normTypes.get(0));
        assertEquals("VC", normTypes.get(1));
        assertEquals("VC_SQRT", normTypes.get(2));
        assertEquals("KR", normTypes.get(3));
    }

    @Test
    public void getNVI() {
        String nvi = NVI.getNVI("https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/e351f7cc-7a2c-4515-ae0b-3bb2f91c986a/4DNFIMIMLMD3.hic");
        assertEquals("240094740,25900", nvi);
    }

}
