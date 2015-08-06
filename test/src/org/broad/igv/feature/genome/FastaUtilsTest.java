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

package org.broad.igv.feature.genome;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertNotNull;

/**
 * User: jacob
 * Date: 2013-Mar-01
 */
public class FastaUtilsTest extends AbstractHeadlessTest {

    @Test
    public void testCreateIndex_01() throws Exception {
        String inPath = TestUtils.DATA_DIR + "fasta/ecoli_out.padded2.fasta";
        String outPath = tstCreateIndex(inPath);

        FastaIndex index = new FastaIndex(outPath);
        assertEquals(1, index.getSequenceNames().size());
        String contig = "NC_000913_bb";
        assertNotNull(index.getIndexEntry(contig));
    }

    public void tstCreateIndex_02(String inPath) throws Exception {
        String outPath = tstCreateIndex(inPath);

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

        //We assume that the test file will have LF line endings
        //regardless of the platform it's on. May change this in the future.
        //git default checkouts may change line endings
        int bytesAtEnd = TestUtils.getBytesAtEnd(inPath);

        int tAbytesPL = tAbasesPL + bytesAtEnd;

        assertEquals(tAbasesPL, entry.getBasesPerLine());
        assertEquals(tAbytesPL, entry.getBytesPerLine());
        assertEquals(9 + bytesAtEnd, entry.getPosition());
        int tAsize = 7 * tAbasesPL + 29;
        assertEquals(tAsize, entry.getSize());
        assertEquals(tA, entry.getContig());

        entry = index.getIndexEntry(tG);

        int tGbasesPL = 56;
        int tGbytesPL = tGbasesPL + bytesAtEnd;
        assertEquals(tGbasesPL, entry.getBasesPerLine());
        assertEquals(tGbytesPL, entry.getBytesPerLine());
        //Starting position is from tA start + tA length + length of header line
        //Since "size" is number of bases, and "position" is bytes, this
        //may look weird
        long tGpos = tAsize + 8*bytesAtEnd + index.getIndexEntry(tA).getPosition() + 9 + bytesAtEnd;
        assertEquals(tGpos, entry.getPosition());
        int tGsize = 5 * tGbasesPL + 26;
        assertEquals(tGsize, entry.getSize());
        assertEquals(tG, entry.getContig());

        GenomeManager manager = GenomeManager.getInstance();
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
        tstCreateIndex(inPath);
    }

    @Test(expected = DataLoadException.class)
    public void testCreateIndexBlankLines() throws Exception {
        String inPath = TestUtils.DATA_DIR + "fasta/blank_lines.fas";
        tstCreateIndex(inPath);
    }

    @Test
    public void testCreateIndexEcoli() throws Exception {
        String inPath = TestUtils.LARGE_DATA_DIR + "ecoli.fasta";
        String outPath = tstCreateIndex(inPath);

        GenomeManager manager = GenomeManager.getInstance();
        Genome genome = manager.loadGenome(inPath, null);
        String chr = "gi|110640213|ref|NC_008253.1|";
        assertNotNull(genome.getChromosome(chr));
        //See http://www.ncbi.nlm.nih.gov/nuccore/110640213
        assertEquals(4938920, genome.getTotalLength());

        String beg = "ATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAG";
        int begloc = 30 - 1;
        int endloc = begloc + beg.length();

        byte[] seq = genome.getSequence(chr, begloc, endloc);
        String sseq = new String(seq);
        assertEquals(beg, sseq);
    }

    static Map<String, Integer> testFastaBlanks;

    static {
        testFastaBlanks = new HashMap<String, Integer>(3);
        testFastaBlanks.put("Chrom1", 632);
        testFastaBlanks.put("Chrom2", 284);
        testFastaBlanks.put("Chrom3", 287);
    }

    @Test
    public void testCreateIndexTrailingBlank() throws Exception {
        String inPath = TestUtils.DATA_DIR + "fasta/trailing_line.fas";
        tstCreateIndexGoodBlanks(inPath, testFastaBlanks);
    }

    @Test
    public void testCreateIndexBlankBetweenContigs() throws Exception {
        String inPath = TestUtils.DATA_DIR + "fasta/blank_lines_betweencontigs.fas";
        tstCreateIndexGoodBlanks(inPath, testFastaBlanks);
    }


    private String tstCreateIndex(String inPath) throws Exception{
        String outPath = inPath + ".fai";
        File outFile = new File(outPath);
        outFile.delete();
        outFile.deleteOnExit();

        FastaUtils.createIndexFile(inPath, outPath);
        return outPath;
    }
    public void tstCreateIndexGoodBlanks(String inPath, Map<String, Integer> expectedSizes) throws Exception {
        String outPath = tstCreateIndex(inPath);

        FastaIndex seq = new FastaIndex(outPath);
        assertEquals(expectedSizes.size(), seq.getSequenceNames().size());
        for (String chrom : expectedSizes.keySet()) {
            assertEquals((int) expectedSizes.get(chrom), seq.getSequenceSize(chrom));
        }
    }

    @Test(expected = DataLoadException.class)
    public void testCreateIndexDuplicateContigs() throws Exception{
        String inPath = TestUtils.DATA_DIR + "fasta/dup_contigs.fas";
        String outPath = TestUtils.DATA_DIR + "out/tmp.fai";
        File outFile = new File(outPath);
        outFile.delete();
        outFile.deleteOnExit();

        FastaUtils.createIndexFile(inPath, outPath);

    }

}
