/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.feature.genome;

import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.junit.*;

import java.io.IOException;

import static org.junit.Assert.assertEquals;

/**
 * @author jrobinso
 */
public class SequenceHelperTest {

    static String seqPath = "http://igvdata.broadinstitute.org/genomes/seq/hg18/";
    static PreferenceManager preferenceManager;
    static boolean useByteRange;

    public SequenceHelperTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        Globals.setHeadless(true);
        preferenceManager = PreferenceManager.getInstance();
        useByteRange = preferenceManager.getAsBoolean(PreferenceManager.USE_BYTE_RANGE);
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
        preferenceManager.put(PreferenceManager.USE_BYTE_RANGE, useByteRange);
    }

    @Before
    public void setUp() {
        //Web requests don't seem to work with this false
        preferenceManager.put(PreferenceManager.USE_BYTE_RANGE, true);
    }

    @After
    public void tearDown() {
        SequenceHelper.setCacheSequences(true);
    }

    /**
     * Test cached vs uncached sequence reads.
     */
    @Test
    public void testCache() {
        try {
            String chr = "chr7";
            int start = 55054464;
            int end = start + 10000;


            SequenceHelper helper = new SequenceHelper(seqPath);

            byte[] cachedSeq = helper.getSequence(chr, start, end, Integer.MAX_VALUE);

            SequenceHelper.setCacheSequences(false);

            byte[] uncachedSeq = helper.getSequence(chr, start, end, Integer.MAX_VALUE);

            assertEquals(uncachedSeq.length, cachedSeq.length);

            for (int i = 0; i < cachedSeq.length; i++) {
                assertEquals("i=" + i, uncachedSeq[i], cachedSeq[i]);
            }

        } finally {
            SequenceHelper.setCacheSequences(true);
        }

    }

    @Test
    public void testNonCloudGenome() throws IOException {


        String chr = "scaffold_v2_26164";
        int start = 5;
        int end = 10;
        String expSequence = "ATTGC";
        SequenceHelper helper = new SequenceHelper("http://www.broadinstitute.org/igvdata/annotations/seq/spur_2.1/");
        byte[] seq = helper.getSequence(chr, start, end, Integer.MAX_VALUE);
        assertEquals(expSequence, new String(seq));

    }

    /**
     * Test known sequence (start of EGFR).
     */
    @Test
    public void readEGFRSequence() {
        String chr = "chr7";
        int start = 55054464;
        int end = start + 20;
        String expSequence = "ATGCGACCCTCCGGGACGGC";
        SequenceHelper helper = new SequenceHelper(seqPath);
        byte[] seq = helper.getSequence(chr, start, end, Integer.MAX_VALUE);
        assertEquals(expSequence, new String(seq));
    }


    @Test
    public void testByteRangeFalse() {
        PreferenceManager preferenceManager = PreferenceManager.getInstance();
        preferenceManager.put(PreferenceManager.USE_BYTE_RANGE, false);

        readEGFRSequence();
    }

}