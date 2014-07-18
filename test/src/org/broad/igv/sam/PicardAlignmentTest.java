/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
* To change this template, choose Tools | Templates
* and open the template in the editor.
 */

package org.broad.igv.sam;

//~--- non-JDK imports --------------------------------------------------------

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
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
