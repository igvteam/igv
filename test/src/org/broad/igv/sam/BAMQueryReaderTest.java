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

package org.broad.igv.sam;

import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.Globals;
import org.broad.igv.sam.reader.BAMQueryReader;
import org.broad.igv.sam.reader.SamQueryTextReader;
import org.broad.igv.util.TestUtils;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * @author jrobinso
 */
public class BAMQueryReaderTest {

    public BAMQueryReaderTest() {
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
    @Test
    public void testClose() throws Exception {
        assertTrue(false);
    }

    /**
     * Test of query method, of class BAMQueryReader.
     */
    @Test
    public void testQueryAgainstSam() throws Exception {

        String bamfile = TestUtils.LARGE_DATA_DIR + "/HG00171.hg18.bam";
        String samfile = TestUtils.LARGE_DATA_DIR + "/HG00171.hg18.sam";
        String chr = "chr1";
        int end = 300000000;
        int start = end / 5;

        BAMQueryReader bamreader = new BAMQueryReader(new File(bamfile));
        SamQueryTextReader samreader = new SamQueryTextReader(samfile);
        CloseableIterator<Alignment> bamiter = bamreader.query(chr, start, end, true);
        CloseableIterator<Alignment> samiter = samreader.iterator();
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

}