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

package org.broad.igv.sam;

import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.sam.reader.*;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author jrobinso
 */
public class SamIndexerTest extends AbstractHeadlessTest{

    /**
     * Test of createSamIndex method, of class SamIndexer.
     */
    @Test
    public void createSamIndex() throws Exception {
        String samFilePath = TestUtils.DATA_DIR + "sam/test_2.sam";
        File samFile = new File(samFilePath);
        File idxFile = new File(samFilePath + ".sai");
        idxFile.deleteOnExit();
        int tileWidth = 100;
        AlignmentIndexer instance = new SamIndexer(samFile, null, null);
        //AlignmentIndexer.getInstance(samFile, null, null);

        FeatureIndex result = instance.createSamIndex(idxFile, tileWidth);
        result.store(idxFile);

        SAMReader reader = new SAMReader(samFilePath, true);

        String chr = "chr3";
        int start = 125963669 - 100;
        int end = start + 300;

        assertAlignmentQueryValid(reader, chr, start, end);
    }

    @Test
    public void createBamIndex() throws Exception{
        String bamFilePath = TestUtils.DATA_DIR + "bam/chr1_chr2.hg18.bam";
        File bamFile = new File(bamFilePath);

        String bamIndexPath = bamFilePath + ".bai";
        File bamIndex = new File(bamIndexPath);
        bamIndex.deleteOnExit();
        SamIndexer.createBAMIndex(bamFile, bamIndex);

        BAMFileReader reader = new BAMFileReader(bamFile);

        assertTrue("No index found", reader.hasIndex());

        String chr = "chr1";
        int start = 155743955 - 100;
        int end = start + 300;

        assertAlignmentQueryValid(reader, chr, start, end);
    }

    private void assertAlignmentQueryValid(AlignmentReader reader, String chr, int start, int end) throws IOException{

        CloseableIterator<Alignment> iter = reader.query(chr, start, end, false);
        int count = 0;

        while (iter.hasNext()) {
            Alignment rec = iter.next();
            assertEquals(chr, rec.getChr());
            assertTrue(rec.getStart() >= start);
            assertTrue(rec.getStart() < end);
            count++;
        }

        iter.close();

        //System.out.println(count + " alignments returned by query");
        assertTrue("No alignments returned by query", count > 0);

    }
}
