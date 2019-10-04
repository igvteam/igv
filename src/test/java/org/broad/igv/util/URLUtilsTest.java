package org.broad.igv.util;

import org.junit.Test;

import java.net.MalformedURLException;

import static org.junit.Assert.*;

public class URLUtilsTest {

    @Test
    public void getPath() throws MalformedURLException {
        String path = URLUtils.getPath("gs://foo.bar/test.bam?a=b");
        assertEquals("/test.bam", path);
    }

    @Test
    public void isURL() {
        assertTrue(URLUtils.isURL("https://foo.bar"));
        assertTrue(URLUtils.isURL("http://foo.bar"));
        assertTrue(URLUtils.isURL("gs://foo.bar"));
        assertTrue(URLUtils.isURL("s3://foo.bar"));
        assertTrue(URLUtils.isURL("file://foo.bar"));
    }

    @Test
    public void addExtension() {
        String extension = ".bai";
        String url = "http://foo/bar.bam?key1=a&key2=b";
        String ammendedURL = "http://foo/bar.bam.bai?key1=a&key2=b";
        assertEquals(ammendedURL, URLUtils.addExtension(url, extension));
    }

    @Test
    public void replaceExtension() {
        String extension = ".bai";
        String url = "http://foo/bar.bam?key1=a&key2=b";
        String ammendedURL = "http://foo/bar.bai?key1=a&key2=b";
        assertEquals(ammendedURL, URLUtils.replaceExtension(url, ".bam", extension));

    }

    @Test
    public void getQuery() throws MalformedURLException {

        String q = URLUtils.getQuery("http://noquery");
        assertNull(q);

        q = URLUtils.getQuery("s3://foo/bar.bam?key1=a&key2=b");
        assertEquals("key1=a&key2=b", q);
    }

}