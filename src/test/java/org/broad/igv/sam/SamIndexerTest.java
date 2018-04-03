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
