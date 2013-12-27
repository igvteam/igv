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

import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.Globals;
import org.broad.igv.sam.reader.BAMFileReader;
import org.broad.igv.sam.reader.SAMReader;
import org.broad.igv.util.TestUtils;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import java.io.File;

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

    @Ignore("Just scratch code")
    @Test
    public void testCopyBAMFileSnippet() {

        String sequence = "chr1";
        int end = 300000000;
        int start = end / 5;

        String bamFile = TestUtils.LARGE_DATA_DIR + "HG00171.hg18.bam";
        SAMFileReader reader = new SAMFileReader(new File(bamFile));
        reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        // Hit the index to determine the chunk boundaries for the required data.
        final SAMFileHeader fileHeader = reader.getFileHeader();
        final int referenceIndex = fileHeader.getSequenceIndex(sequence);

        String outPath = "tmpbam.bam";
        File outFile = new File(outPath);
        File indexFile = new File(outPath.replace(".bam", ".bai"));
        outFile.delete();
        indexFile.delete();
        outFile.deleteOnExit();
        indexFile.deleteOnExit();
        SAMFileWriterFactory factory = new SAMFileWriterFactory();
        //May wish to tinker with options
        //factory.setMaxRecordsInRam(1000).setUseAsyncIo(true);
        boolean createIndex = true;
        factory.setCreateIndex(createIndex);
        SAMFileWriter writer = factory.makeSAMOrBAMWriter(fileHeader, true, outFile);

        int count = 0;
        //NOTE: 1-BASED START/END
        SAMRecordIterator iter = reader.queryOverlapping(sequence, start, end);
        while (iter.hasNext()) {
            writer.addAlignment(iter.next());
            count++;
        }
        iter.close();
        writer.close();


        assertEquals("Index file existence unexpected: " + indexFile.getAbsolutePath(), createIndex, indexFile.exists());

        SAMFileReader writtenReader = new SAMFileReader(new File(outPath));
        writtenReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        iter = null;
        if(createIndex){
            iter = writtenReader.queryOverlapping(sequence, start, end);
        }else{
            iter = writtenReader.iterator();
        }

        int readCount = 0;
        while (iter.hasNext()) {
            readCount++;
            iter.next();
        }

        System.out.println(readCount + " alignments read");
        assertTrue("No alignments read", readCount > 0);
        assertEquals("Read a different number of alignments than written", count, readCount);


    }

}