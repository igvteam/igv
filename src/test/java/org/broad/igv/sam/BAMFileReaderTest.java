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
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.sam.reader.BAMReader;
import org.broad.igv.sam.reader.SAMReader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * @author jrobinso
 */
public class BAMFileReaderTest {

    public BAMFileReaderTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        Globals.setHeadless(true);
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    /**
     * Test of close method, of class BAMQueryReader.
     */
    @Test @Ignore("Requires largedata bundle")
    public void testClose() throws Exception {
        String bamfile = TestUtils.LARGE_DATA_DIR + "HG00171.hg18.bam";
        String chr = "chr1";
        int end = 300000000;
        int start = end / 5;
        int stopafter = 10;
        int counter = 0;
        BAMReader bamreader = new BAMReader(new ResourceLocator(bamfile), true);
        CloseableIterator<SAMAlignment> bamiter = bamreader.query(chr, start, end, true);
        while (bamiter.hasNext()) {
            Alignment bamrecord = bamiter.next();
            if (counter >= stopafter) {
                break;
            } else {
                counter++;
            }
        }
        bamreader.close();
        boolean closeSucceeded = false;
        try {
            CloseableIterator<SAMAlignment> bamiter2 = bamreader.query(chr, start, end, true);
        } catch (NullPointerException npe) {
            closeSucceeded = true;
        }
        assertTrue(closeSucceeded);
    }

    /**
     * Test of query method, of class BAMQueryReader.
     */
    @Test @Ignore("Requires largedata bundle")
    public void testQueryAgainstSam() throws Exception {

        String bamfile = TestUtils.LARGE_DATA_DIR + "HG00171.hg18.bam";
        String samfile = TestUtils.LARGE_DATA_DIR + "HG00171.hg18.sam";
        String chr = "chr1";
        int end = 300000000;
        int start = end / 5;

        BAMReader bamreader = new BAMReader(new ResourceLocator(bamfile), true);
        SAMReader samreader = new SAMReader(samfile);
        CloseableIterator<SAMAlignment> bamiter = bamreader.query(chr, start, end, true);
        CloseableIterator<SAMAlignment> samiter = samreader.iterator();
        int count = 0;
        while (bamiter.hasNext()) {
            Alignment bamrecord = bamiter.next();
            Alignment samrecord = samiter.next();
            assertTrue(bamrecord.getStart() >= start);
            assertTrue(bamrecord.getEnd() <= end);
            assertEquals(bamrecord.getReadName(), samrecord.getReadName());
            assertEquals(bamrecord.getSample(), samrecord.getSample());
            count++;
        }
        assertTrue("No data retrieved", count > 0);
        System.out.println("Retrieved " + count + " rows");

    }

    @Test
    public void testCSI() throws Exception {

        String bamfile = TestUtils.DATA_DIR + "bam/BAMFileIndexTest/index_test.bam";
        String baifile = TestUtils.DATA_DIR + "bam/BAMFileIndexTest/index_test.bam.bai";
        String csifile = TestUtils.DATA_DIR + "bam/BAMFileIndexTest/index_test.bam.csi";

        String chr = "chr1";
        int end = 6000000;
        int start = 1000000;

        ResourceLocator baiLocator = new ResourceLocator(bamfile);
        baiLocator.setIndexPath(baifile);
        BAMReader baireader = new BAMReader(baiLocator, true);

        ResourceLocator csiLocator = new ResourceLocator(bamfile);
        csiLocator.setIndexPath(csifile);
        BAMReader csiReader = new BAMReader(csiLocator, true);



        CloseableIterator<SAMAlignment> baiiter = baireader.query(chr, start, end, true);
        CloseableIterator<SAMAlignment> csiiter = csiReader.query(chr, start, end, true);

        int count = 0;
        while (baiiter.hasNext()) {
            Alignment bamrecord = baiiter.next();
            Alignment samrecord = csiiter.next();
            assertTrue(bamrecord.getStart() >= start);
            assertTrue(bamrecord.getEnd() <= end);
            assertEquals(bamrecord.getReadName(), samrecord.getReadName());
            assertEquals(bamrecord.getSample(), samrecord.getSample());
            count++;
        }
        assertTrue("Unexpected data count: " + count, count == 20);

    }



    public int count(CloseableIterator<SAMAlignment> iter) {
        int counter = 0;
        while (iter.hasNext()) {
            iter.next();
            counter++;
        }
        return counter;
    }

}