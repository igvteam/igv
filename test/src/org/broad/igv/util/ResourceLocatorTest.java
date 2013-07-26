package org.broad.igv.util;


import org.junit.Test;

import java.net.URL;

import static org.junit.Assert.*;


public class ResourceLocatorTest {


    @Test
    public void testGetTypeString() throws Exception {

        // A URL from GenomeBridge
        ResourceLocator loc = new ResourceLocator("https://cmi-jobs-ci.s3.amazonaws.com/75b9e3ad-fafe-439e-8d37-db49acdb7875-1%2F9274485a-20c7-4717-b552-c035a47a03ef.clean.dedup.recal.bam?AWSAccessKeyId=AKIAJ77WFKEZDBZ4ND5Q&Expires=1372349359&Signature=WE1XdJjbQPBWhm3jAFmrBPYVRdA%3D");
        String type = loc.getTypeString();
        assertTrue(type.endsWith(".bam"));

        // Genome space with type converstion
        loc = new ResourceLocator("https://dmtest.genomespace.org:8444/datamanager/files/users/SAGDemo/Step1/TF.data.tab?dataformat=http://www.genomespace.org/datamanager/dataformat/gct/0.0.0");
        type = loc.getTypeString();
        assertTrue(type.endsWith(".gct"));

        // Bam file with non-standard url
        loc = new ResourceLocator("http://some.server.org/foo?file=/server/local/path/bar&param2=value2&dataformat=.bam");
        type = loc.getTypeString();
        assertTrue(type.endsWith(".bam"));

        // Bam file specified by file parameter
        loc = new ResourceLocator("http://some.server.org/foo?file=/server/local/path/bar.bam");
        type = loc.getTypeString();
        assertTrue(type.endsWith(".bam"));


        // Gzipped file
        loc = new ResourceLocator("/foo/bar.cn.gz");
        type = loc.getTypeString();
        assertTrue(type.endsWith(".cn"));

        // File with "double" extension and mixed case
        loc = new ResourceLocator("/foo/bar.PEAK.bin");
        type = loc.getTypeString();
        assertTrue(type.endsWith(".peak.bin"));
    }

    @Test
    public void testBamIndexPaths() throws Exception {

        String url = "http://some.server.org/foo/bar.bam";
        ResourceLocator rl = new ResourceLocator(url);
        String indexPath = rl.getBamIndexPath();
        assertEquals(url + ".bai", indexPath);

        url = "http://some.server.org/foo?file=/server/local/path/bar.bam&param2=value2";
        String expectedPath = "http://some.server.org/foo?file=/server/local/path/bar.bam.bai&param2=value2";
        rl = new ResourceLocator(url);
        indexPath = rl.getBamIndexPath();
        assertEquals(expectedPath, indexPath);

        // Bam with bam extension
        url = "http://some.server.org/foo?file=/server/local/path/bar&param2=value2&dataformat=.bam";
        expectedPath = "http://some.server.org/foo?file=/server/local/path/bar.bai&param2=value2&dataformat=.bam";
        rl = new ResourceLocator(url);
        indexPath = rl.getBamIndexPath();
        assertEquals(expectedPath, indexPath);


        String localPath = "/foo/bar.bam";
        rl = new ResourceLocator(localPath);
        indexPath = rl.getBamIndexPath();
        assertEquals(localPath + ".bai", indexPath);


    }
}
