/*
 *  The MIT License (MIT)
 *
 * Copyright (c) 2016 University of California San Diego
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 *
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.sam.cram;

import htsjdk.samtools.util.CloseableIterator;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.PicardAlignment;
import org.broad.igv.sam.reader.BAMReader;
import org.broad.igv.sam.reader.SAMReader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * @author jrobinso
 */
public class CRAMReaderTest {

    public CRAMReaderTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        Globals.setHeadless(true);
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }


    @Test
    public void testIterateLocalCraiCram() throws Exception {

        String cramFile = TestUtils.DATA_DIR + "cram/cram_with_crai_index.cram";
        GenomeManager.getInstance().loadGenome(TestUtils.DATA_DIR + "cram/hg19mini.fasta", null);

        BAMReader reader = new BAMReader(new ResourceLocator(cramFile), true);

        List<String> seqNames = reader.getSequenceNames();
        assertEquals(4, seqNames.size());


        CloseableIterator<PicardAlignment> iter = reader.iterator();
        int counter = count(iter);
        assertEquals(11, counter);

    }

    @Test
    public void testQueryLocalCraiCram() throws Exception {

        String cramFile = TestUtils.DATA_DIR + "cram/cram_with_crai_index.cram";
        GenomeManager.getInstance().loadGenome(TestUtils.DATA_DIR + "cram/hg19mini.fasta", null);

        BAMReader reader = new BAMReader(new ResourceLocator(cramFile), true);

        CloseableIterator<PicardAlignment> iter = reader.query("2", 500, 600, false);
        int counter = count(iter);
        assertEquals(2, counter);
    }


    @Test
    public void testQueryLocalBaiCram() throws Exception {

        String cramFile = TestUtils.DATA_DIR + "cram/cram_with_bai_index.cram";
        GenomeManager.getInstance().loadGenome(TestUtils.DATA_DIR + "cram/hg19mini.fasta", null);

        BAMReader reader = new BAMReader(new ResourceLocator(cramFile), true);

        CloseableIterator<PicardAlignment> iter = reader.query("2", 500, 600, false);
        int counter = count(iter);
        assertEquals(2, counter);
    }



    @Test
    public void testRemoteCraiCram() throws Exception {

        String cramFile = "https://s3.amazonaws.com/igv.broadinstitute.org/test/cram/cram_with_crai_index.cram";
        GenomeManager.getInstance().loadGenome(TestUtils.DATA_DIR + "cram/hg19mini.fasta", null);

        BAMReader reader = new BAMReader(new ResourceLocator(cramFile), true);

        CloseableIterator<PicardAlignment> iter = reader.query("2", 500, 600, false);
        int counter = count(iter);
        assertEquals(2, counter);
    }


    @Test
    public void testRemoteBaiCram() throws Exception {

        String cramFile = TestUtils.DATA_DIR + "cram/cram_with_bai_index.cram";
        GenomeManager.getInstance().loadGenome(TestUtils.DATA_DIR + "cram/hg19mini.fasta", null);

        BAMReader reader = new BAMReader(new ResourceLocator(cramFile), true);

        CloseableIterator<PicardAlignment> iter = reader.query("2", 500, 600, false);
        int counter = count(iter);
        assertEquals(2, counter);
    }


    public int count(CloseableIterator<PicardAlignment> iter) {
        int counter = 0;
        while (iter.hasNext()) {
            iter.next();
            counter++;
        }
        return counter;
    }

}