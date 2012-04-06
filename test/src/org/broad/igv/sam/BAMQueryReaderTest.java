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

import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.Globals;
import org.broad.igv.sam.reader.BAMFileReader;
import org.broad.igv.sam.reader.SAMReader;
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
        String bamfile = TestUtils.LARGE_DATA_DIR + "HG00171.hg18.bam";
        String chr = "chr1";
        int end = 300000000;
        int start = end / 5;
        int stopafter = 10;
        int counter = 0;
        BAMFileReader bamreader = new BAMFileReader(new File(bamfile));
        CloseableIterator<Alignment> bamiter = bamreader.query(chr, start, end, true);
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
            CloseableIterator<Alignment> bamiter2 = bamreader.query(chr, start, end, true);
        } catch (NullPointerException npe) {
            closeSucceeded = true;
        }
        assertTrue(closeSucceeded);


    }

    /**
     * Test of query method, of class BAMQueryReader.
     */
    @Test
    public void testQueryAgainstSam() throws Exception {

        String bamfile = TestUtils.LARGE_DATA_DIR + "HG00171.hg18.bam";
        String samfile = TestUtils.LARGE_DATA_DIR + "HG00171.hg18.sam";
        String chr = "chr1";
        int end = 300000000;
        int start = end / 5;

        BAMFileReader bamreader = new BAMFileReader(new File(bamfile));
        SAMReader samreader = new SAMReader(samfile);
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