package org.broad.igv.sam.reader;

import net.sf.samtools.util.CloseableIterator;
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
