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

package org.broad.igv.feature.genome;

import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.util.TestUtils;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import static junit.framework.Assert.*;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: 8/7/11
 * Time: 7:41 PM
 * To change this template use File | Settings | File Templates.
 */
public class FastaIndexTest {

    static String indexPath = "http://www.broadinstitute.org/igvdata/test/fasta/ci2_test.fa.fai";
    private FastaIndex index;

    @Before
    public void setUp() throws IOException {
        TestUtils.setUpHeadless();

        index = new FastaIndex(indexPath);
    }

    @Test
    public void testGetContigs() throws Exception {

        List expectedContigs = Arrays.asList("chr01p", "chr02q", "chr03q");

        Set<String> contigs = index.getSequenceNames();
        assertEquals(expectedContigs.size(), contigs.size());
        for (Object exp : expectedContigs) {
            assertTrue(contigs.contains(exp));
        }

    }

    @Test
    public void testGetIndexEntry() throws Exception {
        //5803340	14565324	50	51
        FastaIndex.FastaSequenceIndexEntry entry = index.getIndexEntry("chr03q");
        assertEquals("Size", 5803340, entry.getSize());
        assertEquals("Position", 14565324, entry.getPosition());
        assertEquals("basesPerLine", 50, entry.getBasesPerLine());
        assertEquals("bytesPerLine", 51, entry.getBytesPerLine());
    }


    @Test
    public void testCreateIndex_01() throws Exception {
        String inPath = TestUtils.DATA_DIR + "fasta/ecoli_out.padded2.fasta";
        String outPath = TestUtils.DATA_DIR + "out/ecoli_out.padded2.fasta.fai";
        File outFile = new File(outPath);
        outFile.delete();
        outFile.deleteOnExit();

        FastaIndex.createIndexFile(inPath, outPath);

        FastaIndex index = new FastaIndex(outPath);
        assertEquals(1, index.getSequenceNames().size());
        String contig = "NC_000913_bb";
        assertNotNull(index.getIndexEntry(contig));
    }

    public void tstCreateIndex_02(String inPath) throws Exception {
        String outPath = inPath + ".fai";
        File outFile = new File(outPath);
        outFile.delete();
        outFile.deleteOnExit();

        FastaIndex.createIndexFile(inPath, outPath);

        FastaIndex index = new FastaIndex(outPath);
        assertEquals(2, index.getSequenceNames().size());
        String tA = "my:testA";
        String tG = "my:testG";
        String[] contigs = {tA, tG};
        for (String contig : contigs) {
            assertNotNull(index.getIndexEntry(contig));
        }

        FastaIndex.FastaSequenceIndexEntry entry = index.getIndexEntry(tA);
        int tAbasesPL = 58;

        int tAbytesPL = tAbasesPL + 1;

        assertEquals(tAbasesPL, entry.getBasesPerLine());
        assertEquals(tAbytesPL, entry.getBytesPerLine());
        assertEquals(10, entry.getPosition());
        int tAsize = 7 * tAbasesPL + 29;
        assertEquals(tAsize, entry.getSize());
        assertEquals(tA, entry.getContig());

        entry = index.getIndexEntry(tG);

        int tGbasesPL = 56;
        int tGbytesPL = tGbasesPL + 1;
        assertEquals(tGbasesPL, entry.getBasesPerLine());
        assertEquals(tGbytesPL, entry.getBytesPerLine());
        //Starting position is from tA start + tA length + length of header line
        //Since "size" is number of bases, and "position" is bytes, this
        //may look weird
        long tGpos = tAsize + 8 + index.getIndexEntry(tA).getPosition() + 10;
        assertEquals(tGpos, entry.getPosition());
        int tGsize = 5 * tGbasesPL + 26;
        assertEquals(tGsize, entry.getSize());
        assertEquals(tG, entry.getContig());

        GenomeManager manager = new GenomeManager();
        Genome genome = manager.loadGenome(inPath, null);
        String tAseq = new String(genome.getSequence(tA, 0, tAsize));
        assertEquals(tAsize, tAseq.length());
        String remmed = tAseq.replaceAll("A", "");
        assertEquals(0, remmed.length());

        String tGseq = new String(genome.getSequence(tG, 0, tGsize));
        assertEquals(tGsize, tGseq.length());
        remmed = tGseq.replaceAll("G", "");
        assertEquals(0, remmed.length());
    }


    @Test
    public void testCreateIndexUncompressed() throws Exception {
        String inPath = TestUtils.DATA_DIR + "fasta/fasta_2contigs.fa";
        tstCreateIndex_02(inPath);
    }

    @Test(expected = DataLoadException.class)
    public void testCreateIndexUneven() throws Exception {
        String inPath = TestUtils.DATA_DIR + "fasta/fasta_uneven.fa";
        String outPath = TestUtils.DATA_DIR + "out/tmp.fai";
        File outFile = new File(outPath);
        outFile.delete();
        outFile.deleteOnExit();

        FastaIndex.createIndexFile(inPath, outPath);
    }

    @Test
    public void testCreateIndexEcoli() throws Exception {
        String inPath = TestUtils.LARGE_DATA_DIR + "ecoli.fasta";
        String outPath = inPath + ".fai";
        File outFile = new File(outPath);

        outFile.delete();
        outFile.deleteOnExit();

        GenomeManager manager = new GenomeManager();
        Genome genome = manager.loadGenome(inPath, null);
        String chr = "gi|110640213|ref|NC_008253.1|";
        assertNotNull(genome.getChromosome(chr));
        //See http://www.ncbi.nlm.nih.gov/nuccore/110640213
        assertEquals(4938920, genome.getLength());

        String beg = "ATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAG";
        int begloc = 30 - 1;
        int endloc = begloc + beg.length();

        byte[] seq = genome.getSequence(chr, begloc, endloc);
        String sseq = new String(seq);
        assertEquals(beg, sseq);
    }

}
