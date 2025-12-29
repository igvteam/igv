package org.igv.util;

import org.junit.Test;

import static org.junit.Assert.*;

public class GoogleUtilsTest {

    @Test
    public void translateGoogleCloudURL() {
        String gsURL = "gs://genomics-public-data/platinum-genomes/bam/NA12877_S1.bam";
        String url = GoogleUtils.translateGoogleCloudURL(gsURL);
        assertEquals("https://storage.googleapis.com/storage/v1/b/genomics-public-data/o/platinum-genomes%2Fbam%2FNA12877_S1.bam?alt=media", url);

       assertTrue(GoogleUtils.isGoogleStorageURL(url));

    }

    @Test
    public void testEncoding() {
        String str = " !\"#$&'()*+,/:;=?@[]";
        String encoded = "%20%21\"%23%24%26%27%28%29%2A%2B%2C%2F%3A%3B%3D%3F%40%5B%5D";
        assertEquals(encoded, GoogleUtils.googleObjectEncode(str));
    }


}