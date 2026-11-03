package org.igv.htsget;

import org.junit.Test;

import static org.junit.Assert.*;


public class HtsgetReaderTest {

    @Test
    public void testReadHeader() throws Exception {

        String url = "https://htsget.ga4gh-demo.org/variants/spec-v4.3";
        HtsgetUtils.Metadata metadata = HtsgetUtils.getMetadata(url);
        HtsgetReader reader = HtsgetReader.getReader(metadata);
        byte [] headerBytes  = reader.readHeader();
        assertNotNull(headerBytes);
    }

    @Test
    public void testReadData() throws Exception {

        String url = "https://htsget.ga4gh-demo.org/variants/spec-v4.3";
        String chr = "20";
        int start = 0;
        int end = 10000000;

        HtsgetUtils.Metadata metadata =  HtsgetUtils.getMetadata(url);
        HtsgetReader reader = HtsgetReader.getReader(metadata);


        byte[] data = reader.readData(chr, start, end);
        assertNotNull(data);

    }

    /**
     * BGZF EOF block as a base64 data URL, as returned by some servers for the final block.
     */
    @Test
    public void testDecodeDataURL() {
        String dataURL = "data:application/vnd.ga4gh.bam;base64,H4sIBAAAAAAA/wYAQkMCABsAAwAAAAAAAAAAAA==";
        byte[] bytes = HtsgetReader.decodeDataURL(dataURL);
        assertEquals(28, bytes.length);
        assertEquals((byte) 0x1f, bytes[0]);
        assertEquals((byte) 0x8b, bytes[1]);
    }

}