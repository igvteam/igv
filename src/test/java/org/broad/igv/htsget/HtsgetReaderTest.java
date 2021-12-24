package org.broad.igv.htsget;

import org.junit.Ignore;
import org.junit.Test;

import static org.junit.Assert.*;

public class HtsgetReaderTest {

    @Test
    public void testReadHeader() throws Exception {

        String url = "https://htsget.ga4gh.org/variants/giab.NA12878";
        HtsgetUtils.Metadata metadata = HtsgetUtils.getMetadata(url);
        HtsgetReader reader = HtsgetReader.getReader(metadata);
        byte [] headerBytes  = reader.readHeader();
        assertNotNull(headerBytes);
    }

    @Ignore   // Problems with reference server?
    @Test
    public void testReadData() throws Exception {

        String url = "https://htsget.ga4gh.org/variants/giab.NA12878";
        String chr = "8";
        int start = 128732400 - 1;
        int end = 128770475;

        HtsgetUtils.Metadata metadata =  HtsgetUtils.getMetadata(url);
        HtsgetReader reader = HtsgetReader.getReader(metadata);
        byte[] data = reader.readData(chr, start, end);
        assertNotNull(data);

    }

}