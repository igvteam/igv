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

import org.broad.igv.util.TestUtils;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

/**
 * @author jrobinso
 */
public class SamAlignmentTest {

    /**
     * Constructs ...
     */
    public SamAlignmentTest() {
    }

    /**
     * Method description
     *
     * @throws Exception
     */
    @BeforeClass
    public static void setUpClass() throws Exception {
        TestUtils.setUpHeadless();
    }

    /**
     * Method description
     *
     * @throws Exception
     */
    @AfterClass
    public static void tearDownClass() throws Exception {
    }


    /**
     * Test of adjustReads method, of class SamAlignment.
     */

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

//        SamAlignment instance = new SamAlignment();
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

//        SamAlignment instance = new SamAlignment();
//        instance.createAlignmentBlocks(cigarString, readBases, readBaseQualities);
//
//        byte[] adjustedBases = new byte[80];
//        System.arraycopy(readBases, 0, adjustedBases, 0, 53);
//        adjustedBases[53] = SamAlignment.DELETE_CHAR;
//        adjustedBases[54] = SamAlignment.DELETE_CHAR;
//        adjustedBases[55] = SamAlignment.DELETE_CHAR;
//        adjustedBases[56] = SamAlignment.DELETE_CHAR;
//        System.arraycopy(readBases, 53, adjustedBases, 57, 23);
//
//        assertEquals(adjustedBases.length, instance.getReadLength());
//        for (int i = 0; i < adjustedBases.length; i++)
//        {
//            assertEquals("idx=" + i, adjustedBases[i], instance.getBase(i));
//        }

    }

}
