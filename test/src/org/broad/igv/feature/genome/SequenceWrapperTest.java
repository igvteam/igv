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