/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.sam;

import htsjdk.samtools.util.CloseableIterator;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.sam.reader.*;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author jrobinso
 */
public class SamIndexerTest extends AbstractHeadlessTest{

    /**
     * Test of createSamIndex method, of class SamIndexer.
     */
    @Test
    public void createSamIndex() throws Exception {
        String samFilePath = TestUtils.DATA_DIR + "sam/test_2.sam";
        File samFile = new File(samFilePath);
        File idxFile = new File(samFilePath + ".sai");
        idxFile.deleteOnExit();
        int tileWidth = 100;
        AlignmentIndexer instance = new SamIndexer(samFile, null, null);
        //AlignmentIndexer.getInstance(samFile, null, null);

        FeatureIndex result = instance.createSamIndex(idxFile, tileWidth);
        result.store(idxFile);

        SAMReader reader = new SAMReader(samFilePath, true);

        String chr = "chr3";
        int start = 125963669 - 100;
        int end = start + 300;

        assertAlignmentQueryValid(reader, chr, start, end);
    }



    private void assertAlignmentQueryValid(AlignmentReader reader, String chr, int start, int end) throws IOException{

        CloseableIterator<Alignment> iter = reader.query(chr, start, end, false);
        int count = 0;

        while (iter.hasNext()) {
            Alignment rec = iter.next();
            assertEquals(chr, rec.getChr());
            assertTrue(rec.getStart() >= start);
            assertTrue(rec.getStart() < end);
            count++;
        }

        iter.close();

        //System.out.println(count + " alignments returned by query");
        assertTrue("No alignments returned by query", count > 0);

    }
}
