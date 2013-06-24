package org.broad.igv.util;


import static org.junit.Assert.*;

/**
 * Created with IntelliJ IDEA.
 * User: jrobinso
 * Date: 6/24/13
 * Time: 12:58 PM
 * To change this template use File | Settings | File Templates.
 */
public class ResourceLocatorTest {


    @org.junit.Test
    public void testGetTypeString() throws Exception {

        // A URL from GenomeBridge
        ResourceLocator loc = new ResourceLocator("https://cmi-jobs-ci.s3.amazonaws.com/75b9e3ad-fafe-439e-8d37-db49acdb7875-1%2F9274485a-20c7-4717-b552-c035a47a03ef.clean.dedup.recal.bam?AWSAccessKeyId=AKIAJ77WFKEZDBZ4ND5Q&Expires=1372349359&Signature=WE1XdJjbQPBWhm3jAFmrBPYVRdA%3D");
        String type = loc.getTypeString();
        assertTrue(type.endsWith(".bam"));

        // Genome space with type converstion
        loc = new ResourceLocator("https://dmtest.genomespace.org:8444/datamanager/files/users/SAGDemo/Step1/TF.data.tab?dataformat=http://www.genomespace.org/datamanager/dataformat/gct/0.0.0");
        type = loc.getTypeString();
        assertTrue(type.endsWith(".gct"));

        // Gzipped file
        loc = new ResourceLocator("/foo/bar.cn.gz");
        type = loc.getTypeString();
        assertTrue(type.endsWith(".cn"));
    }
}
