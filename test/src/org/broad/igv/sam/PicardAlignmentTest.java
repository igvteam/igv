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

//~--- non-JDK imports --------------------------------------------------------

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Iterator;

import static junit.framework.Assert.assertEquals;

/**
 * @author jrobinso
 */
public class PicardAlignmentTest extends AbstractHeadlessTest {


    @Ignore
    @Test
    public void testInsertion() {
        byte[] readBases =
                "AATTATTAAAGTGCACTGATCTTTCTTTTCTTTTTCCTTACTATACTTTTTTTTTTTTTTTGAGATGGGAGTTTGG"
                        .getBytes();
        byte[] readBaseQualities =
                "/...$)//2/&/./(/(//+*112,.3/3,11332/113)/6-4446;;;;;;;;;;;;;;;;;;;;;;;;;;;;;"
                        .getBytes();
        String cigarString = "46M1I29M";

//        PicardAlignment instance = new PicardAlignment();
//        instance.createAlignmentBlocks(cigarString, readBases, readBaseQualities);
//
//        byte[] adjustedBases = new byte[75];
//        System.arraycopy(readBases, 0, adjustedBases, 0, 46);
//        System.arraycopy(readBases, 47, adjustedBases, 46, 29);

        /*assertEquals(adjustedBases.length, instance.getBases().length);
        for (int i = 0; i < adjustedBases.length; i++)
        {
            assertEquals("idx=" + i, adjustedBases[i], instance.getBases()[i]);
        }*/

    }

    /**
     * Method description
     */

    @Ignore
    @Test
    public void testDeletion() {
        byte[] readBases =
                "AAGACTCGTGATACTAGCAGAAAATATCAAGAATTAGAATATATATATATATATGTGTGTGTGTGTGTATGTATAC"
                        .getBytes();
        byte[] readBaseQualities =
                "/2231/$13$1,///,/1,/3)31/3/11,/3&*/1$3//367676;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;"
                        .getBytes();
        String cigarString = "53M4D23M";

//        PicardAlignment instance = new PicardAlignment();
//        instance.createAlignmentBlocks(cigarString, readBases, readBaseQualities);
//
//        byte[] adjustedBases = new byte[80];
//        System.arraycopy(readBases, 0, adjustedBases, 0, 53);
//        adjustedBases[53] = PicardAlignment.DELETE_CHAR;
//        adjustedBases[54] = PicardAlignment.DELETE_CHAR;
//        adjustedBases[55] = PicardAlignment.DELETE_CHAR;
//        adjustedBases[56] = PicardAlignment.DELETE_CHAR;
//        System.arraycopy(readBases, 53, adjustedBases, 57, 23);
//
//        assertEquals(adjustedBases.length, instance.getReadLength());
//        for (int i = 0; i < adjustedBases.length; i++)
//        {
//            assertEquals("idx=" + i, adjustedBases[i], instance.getBase(i));
//        }

    }

    /**
     * Test handling of IonTorrent flow signal tags
     * Should be read if available, not cause error if wrong format
     */
    @Test
    public void testFlowSignalTags() throws Exception {
        String inpath = TestUtils.DATA_DIR + "sam/zf_tags.sam";
        AlignmentReader reader = AlignmentReaderFactory.getReader(new ResourceLocator(inpath));
        Iterator<Alignment> iter = reader.iterator();

        int[] expFlowStarts = new int[]{-1, 5, -1};
        int row = 0;
        while (iter.hasNext()) {
            PicardAlignment alignment = (PicardAlignment) iter.next();

            assertEquals(expFlowStarts[row], alignment.getFlowSignalsStart());

            row++;
        }
    }



    @Ignore("Don't store hard clipped bases in blocks")
    @Test
    public void testBuildReadSequenceFromBlocks() throws Exception{
        AlignmentDataManager dataManager = AlignmentDataManagerTest.getManager171();
        Iterator<Alignment> iter = dataManager.getReader().iterator();

        while(iter.hasNext()){
            PicardAlignment al = (PicardAlignment) iter.next();
            assertEquals(al.getRecord().getReadString(), al.buildReadSequenceFromBlocks());
        }
    }
}
