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
public class Uncalled4Test {

    public Uncalled4Test() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        Globals.setHeadless(true);
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Test
    public void testCSI() throws Exception {

        String bamfile = TestUtils.DATA_DIR + "bam/Uncalled_4/dna_r9_dm_ctl.align.bam";
        String baifile = TestUtils.DATA_DIR + "bam/Uncalled_4/dna_r9_dm_ctl.align.bam.bai";


        String chr = "CM044676.1";
        int end = 0;
        int start = Integer.MAX_VALUE;

        ResourceLocator baiLocator = new ResourceLocator(bamfile);
        baiLocator.setIndexPath(baifile);
        BAMReader baireader = new BAMReader(baiLocator, true);

        CloseableIterator<SAMAlignment> baiiter = baireader.query(chr, start, end, true);

        List<String> comments = baireader.getFileHeader().getComments();

        int count = 0;
        while (baiiter.hasNext()) {
            Alignment bamrecord = baiiter.next();

            Object uc = bamrecord.getAttribute("uc");
            Object ud = bamrecord.getAttribute("ud");
            Object ur = bamrecord.getAttribute("ur");
            assertTrue(bamrecord.getStart() >= start);
            assertTrue(bamrecord.getEnd() <= end);
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