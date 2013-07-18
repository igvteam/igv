/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
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

package org.broad.igv.feature.genome;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.PreferenceManager;
import org.broad.igv.util.HttpUtils;
import org.junit.After;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

/**
 * @author jrobinso
 */
public class SequenceWrapperTest extends AbstractHeadlessTest{

    static String seqPath = "http://igvdata.broadinstitute.org/genomes/seq/hg18/";
    static PreferenceManager preferenceManager;
    static SequenceWrapper helper;

    public SequenceWrapperTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        AbstractHeadlessTest.setUpClass();

        preferenceManager = PreferenceManager.getInstance();
        String tmp = SequenceWrapper.checkSequenceURL(seqPath);
        helper = new SequenceWrapper(new IGVSequence(tmp));
    }

    @After
    public void tearDown() throws Exception{
        super.tearDown();
        SequenceWrapper.setCacheSequences(true);
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

            byte[] cachedSeq = helper.getSequence(chr, start, end);

            SequenceWrapper.setCacheSequences(false);

            byte[] uncachedSeq = helper.getSequence(chr, start, end);

            assertEquals(uncachedSeq.length, cachedSeq.length);

            for (int i = 0; i < cachedSeq.length; i++) {
                assertEquals("i=" + i, uncachedSeq[i], cachedSeq[i]);
            }

        } finally {
            SequenceWrapper.setCacheSequences(true);
        }

    }

    @Test
    public void testNonCloudGenome() throws IOException {


        String chr = "scaffold_v2_26164";
        int start = 5;
        int end = 10;
        String expSequence = "ATTGC";
        String tmp = SequenceWrapper.checkSequenceURL("http://www.broadinstitute.org/igvdata/annotations/seq/spur_2.1/");
        SequenceWrapper helper = new SequenceWrapper(new IGVSequence(tmp));
        byte[] seq = helper.getSequence(chr, start, end);
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
        byte[] seq = helper.getSequence(chr, start, end);
        assertEquals(expSequence, new String(seq));
    }


    /**
     * Test fetching sequence using the "Range service" for 3 URL patterns: cloudfront, s3 storage, and the
     * Broad server.
     */
    @Test
    public void testRangeService() {

        try {
            HttpUtils.disableByteRange(true);

            seqPath = "http://igvdata.broadinstitute.org/genomes/seq/hg18/";
            readEGFRSequence();

            seqPath = "http://igv.broadinstitute.org/genomes/seq/hg18/";
            readEGFRSequence();

            seqPath = "http://www.broadinstitute.org/igvdata/annotations/seq/hg18/";
            readEGFRSequence();
        } finally {
            HttpUtils.disableByteRange(false);
        }
    }

    @Test
    public void testGetKeyUniqueness() throws Exception{
        int count = 0;
        int maxNumTiles = 100;
        Set<String> keys = new HashSet<String>();

        for(String chr: genome.getAllChromosomeNames()){
            for(int tNo=0; tNo < maxNumTiles; tNo++){
                String key = SequenceWrapper.getKey(chr, tNo);
                assertFalse("Key not unique: " + key, keys.contains(key));
                keys.add(key);
                count++;
            }
        }

        assertEquals(count, keys.size());

    }

}