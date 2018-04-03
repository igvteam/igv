/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.util;


import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;


public class ResourceLocatorTest {


    @Test
    public void testGetTypeString() throws Exception {

        // A URL from GenomeBridge
        ResourceLocator loc = new ResourceLocator("https://cmi-jobs-ci.s3.amazonaws.com/75b9e3ad-fafe-439e-8d37-db49acdb7875-1%2F9274485a-20c7-4717-b552-c035a47a03ef.clean.dedup.recal.bam?AWSAccessKeyId=AKIAJ77WFKEZDBZ4ND5Q&Expires=1372349359&Signature=WE1XdJjbQPBWhm3jAFmrBPYVRdA%3D");
        String type = loc.getTypeString();
        assertTrue(type.endsWith(".bam"));

        // Genome space with type conversion
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
