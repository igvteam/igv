package org.igv.util;

import org.junit.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;


public class ResourceLocatorTest {


    @Test
    public void testGetFormat() throws Exception {

        // A URL from GenomeBridge
        ResourceLocator loc = new ResourceLocator("https://cmi-jobs-ci.s3.amazonaws.com/75b9e3ad-fafe-439e-8d37-db49acdb7875-1%2F9274485a-20c7-4717-b552-c035a47a03ef.clean.dedup.recal.bam?AWSAccessKeyId=AKIAJ77WFKEZDBZ4ND5Q&Expires=1372349359&Signature=WE1XdJjbQPBWhm3jAFmrBPYVRdA%3D");
        String type = loc.getFormat();
        System.out.println(type);
        assertTrue(type.endsWith("bam"));

        // Bam file with non-standard url
        loc = new ResourceLocator("http://some.server.org/foo?file=/server/local/path/bar&param2=value2&dataformat=bam");
        type = loc.getFormat();
        assertTrue(type.endsWith("bam"));

        // Bam file specified by file parameter
        loc = new ResourceLocator("http://some.server.org/foo?file=/server/local/path/bar.bam");
        type = loc.getFormat();
        assertTrue(type.endsWith("bam"));

        // Gzipped file
        loc = new ResourceLocator("/foo/bar.cn.gz");
        type = loc.getFormat();
        assertTrue(type.endsWith("cn"));


    }

    @Test
    public void getTrackName() throws Exception {

        ResourceLocator resourceLocator = new ResourceLocator("http://foo.bar/track");
        String trackName = resourceLocator.getTrackName();
        assertEquals("track", trackName);

        resourceLocator = new ResourceLocator("http://foo.bar/track?param=abc");
        trackName = resourceLocator.getTrackName();
        assertEquals("track", trackName);

        resourceLocator = new ResourceLocator("http://foo.bar/track?param=abc/123");
        trackName = resourceLocator.getTrackName();
        assertEquals("track", trackName);

    }

    @Test
    public void getLocatorsPairsIndexes() {

        List<File> files = Arrays.asList(
                new File("/data/a.bam"), new File("/data/a.bai"),
                new File("/data/b.cram"), new File("/data/b.crai"),
                new File("/data/c.vcf.gz"), new File("/data/c.vcf.gz.csi"),
                new File("/data/d.vcf"), new File("/data/d.vcf.idx"));

        List<ResourceLocator> locators = ResourceLocator.getLocators(files);

        assertEquals(4, locators.size());
        assertEquals(new File("/data/a.bai").getAbsolutePath(), locators.get(0).getIndexPath());
        assertEquals(new File("/data/b.crai").getAbsolutePath(), locators.get(1).getIndexPath());
        assertEquals(new File("/data/c.vcf.gz.csi").getAbsolutePath(), locators.get(2).getIndexPath());
        assertEquals(new File("/data/d.vcf.idx").getAbsolutePath(), locators.get(3).getIndexPath());
    }

}
