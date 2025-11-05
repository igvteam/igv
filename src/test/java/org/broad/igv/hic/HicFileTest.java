package org.broad.igv.hic;

import htsjdk.samtools.seekablestream.SeekableStream;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.assertEquals;

public class HicFileTest {

    static String testURL = "https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/cc0accc9-a89f-4b40-bac6-de754a7f77c9/4DNFILI6XOB5.hic";
    static HicFile hicFile;

    @BeforeClass
    public static void setUp() throws Exception {
        Map<String, Object> params = new java.util.HashMap<>();
        hicFile =  HicFile.create(testURL, params);
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
        String norm = "NONE";
        List<ContactRecord> records = hicFile.getContactRecords(norm, region1, region2, "BP", binSize, true);
        assertEquals(396963, records.size());
    }

    @Test
    public void getBlocks() {
    }

    @Test
    public void readBlock() {
    }

    @Test
    public void hasNormalizationVector() {
    }

    @Test
    public void getNormalizationVector() {
    }

    @Test
    public void getNormVectorIndex() {
    }

    @Test
    public void readNormVectorIndex() {
    }

    @Test
    public void readNormExpectedValuesAndNormVectorIndex() {
    }

    @Test
    public void skipExpectedValues() {
    }

    @Test
    public void getNormalizationVectorKey() {
    }
}