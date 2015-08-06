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

package org.broad.igv.sam.reader;

import htsjdk.samtools.util.CloseableIterator;
import org.broad.igv.Globals;
import org.broad.igv.sam.Alignment;
import org.junit.BeforeClass;
import org.junit.Test;

import java.net.MalformedURLException;
import java.util.Iterator;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author Jim Robinson
 * @date 1/10/12
 */
public class CGIAlignmentReaderTest {

    @BeforeClass
    public static void setup() {
        Globals.setHeadless(true);
    }

    @Test
    public void testURLs() throws MalformedURLException {

        String testURL = "http://somehost/cramtools/cgi-bin/query.cgi?program=downsample_bam.pl&downsampled_coverage=10&file=input.sam/";
        String queryPath = "query.cgi?";
        String headerPath = "samHeader.cgi?";
        String seqNamePath = "getSequenceNames.cgi?";

        CGIAlignmentReader reader = new CGIAlignmentReader(testURL);

        assertEquals(testURL, reader.getQueryURL());
        assertEquals(testURL.replace(queryPath, headerPath), reader.getHeaderURL());
        assertEquals(testURL .replace(queryPath, seqNamePath), reader.getSequenceNamesURL());

    }

    @Test
     public void testURLsWithPort() throws MalformedURLException {

         String testURL = "http://somehost:8080/cramtools/cgi-bin/query.cgi?program=downsample_bam.pl&downsampled_coverage=10&file=input.sam/";
         String queryPath = "query.cgi?";
         String headerPath = "samHeader.cgi?";
         String seqNamePath = "getSequenceNames.cgi?";

         CGIAlignmentReader reader = new CGIAlignmentReader(testURL);

         assertEquals(testURL, reader.getQueryURL());
         assertEquals(testURL.replace(queryPath, headerPath), reader.getHeaderURL());
         assertEquals(testURL .replace(queryPath, seqNamePath), reader.getSequenceNamesURL());

     }



    // The following test is disabled until we can write our own CGI scripts to test against.
    //@Test
//    public void testQuery() throws Exception {
//
//        String testURL = "http://somehost/query.cgi?file=input.sam";
//
//        CGIAlignmentReader reader = new CGIAlignmentReader(testURL);
//
//        // Assert that we get at least 1 record back.
//        CloseableIterator<Alignment> iter = reader.iterator();
//        assertTrue(iter.hasNext());
//        iter.close();
//
//
//    }


}
