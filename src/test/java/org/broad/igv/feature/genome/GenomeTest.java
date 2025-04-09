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
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.fasta.FastaIndex;
import org.broad.igv.feature.genome.load.GenomeConfig;
import org.broad.igv.util.TestUtils;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import static junit.framework.Assert.*;

/**
 * @author Jim Robinson
 * @date 10/31/11
 */
public class GenomeTest extends AbstractHeadlessTest {

    @Rule
    public TestRule testTimeout = new Timeout((int) 60e3);

    @Test
    public void testGetNCBIName() throws Exception {
        String ncbiID = "gi|125745044|ref|NC_002229.3|";
        String ncbiName = "NC_002229.3";
        assertEquals(ncbiName, Genome.getNCBIName(ncbiID));
    }

    @Test
    public void testGetLongChromosomeNames() throws Exception {
        String mockIndexPath = TestUtils.DATA_DIR + "fasta/bosTau9.fa.fai";
        Sequence sequence = new MockSequence(mockIndexPath);
        GenomeConfig genomeConfig = new GenomeConfig();
        genomeConfig.id = ("bosTau9");
        genomeConfig.setName("bosTau9");
        genomeConfig.setSequence(sequence);
        Genome genome = new Genome(genomeConfig);
        List<String> longChrs = genome.getLongChromosomeNames();
        assertEquals(30, longChrs.size());
    }

    @Test
    public void testGetLongChromosomeNames2() throws Exception {
        String mockIndexPath = TestUtils.DATA_DIR + "fasta/hg19.fa.fai";
        Sequence sequence = new MockSequence(mockIndexPath);
        GenomeConfig genomeConfig = new GenomeConfig();
        genomeConfig.id = ("hg19");
        genomeConfig.setName("hg19");
        genomeConfig.setSequence(sequence);
        Genome genome = new Genome(genomeConfig);
        List<String> longChrs = genome.getLongChromosomeNames();
        assertEquals(24, longChrs.size());
    }

    @Test
    public void testGetLongChromosomeNames3() throws Exception {
        String mockIndexPath = TestUtils.DATA_DIR + "fasta/musa_pseudochromosome.fa.fai";
        Sequence sequence = new MockSequence(mockIndexPath);
        GenomeConfig genomeConfig = new GenomeConfig();
        genomeConfig.id = ("musa_pseudochromosome");
        genomeConfig.setName("musa_pseudochromosome");
        genomeConfig.setSequence(sequence);
        Genome genome = new Genome(genomeConfig);
        List<String> longChrs = genome.getLongChromosomeNames();
        assertEquals(12, longChrs.size());
    }

    @Test
    public void testGetLongChromosomeNames_manySmall() throws Exception {
        String mockIndexPath = TestUtils.DATA_DIR + "fasta/mock_many_small.fa.fai";
        Sequence sequence = new MockSequence(mockIndexPath);
        GenomeConfig genomeConfig = new GenomeConfig();
        genomeConfig.id = ("mock_many_small");
        genomeConfig.setName("mock_many_small");
        genomeConfig.setSequence(sequence);
        Genome genome = new Genome(genomeConfig);
        assertNotNull(genome.getLongChromosomeNames());
        assertTrue("No 'Long' chromosome names found", genome.getLongChromosomeNames().size() > 0);
    }

    /**
     * Class which loads FastaIndex and returns information contained therein,
     * but doesn't actually load full fasta file. For testing
     */
    private class MockSequence implements Sequence {

        private final FastaIndex index;
        private final ArrayList<String> chromoNames;

        public MockSequence(String fastaIndexPath) throws IOException {
            this.index = new FastaIndex(fastaIndexPath);
            this.chromoNames = new ArrayList<String>(index.getSequenceNames());
        }

        @Override
        public byte[] getSequence(String chr, int start, int end) {
            return new byte[0];
        }

        @Override
        public byte getBase(String chr, int position) {
            return 0;
        }

        @Override
        public List<String> getChromosomeNames() {
            return chromoNames;
        }

        @Override
        public int getChromosomeLength(String chrname) {
            return index.getSequenceSize(chrname);
        }

        @Override
        public List<Chromosome> getChromosomes() {
            return index.getChromosomes();
        }

    }

    public static void generateJunkIndex() throws Exception {
        //Generate index file with many small contigs
        int numContigs = 10000;
        int contigMeanSize = 3000;
        int contigSizeRange = 400;
        PrintWriter writer = new PrintWriter(new FileWriter(TestUtils.DATA_DIR + "fasta/mock_many_small.fa.fai"));

        int position = -1;
        int basesPerLine = 80;
        int bytesPerLine = 81;
        for (int ci = 0; ci < numContigs; ci++) {
            String chr = "" + ci;
            int size = contigMeanSize + (int) (contigSizeRange * (Math.random() - 0.5));

            String line = String.format("%s\t%d\t%d\t%d\t%d", chr, size, position, basesPerLine, bytesPerLine);
            writer.println(line);
        }
        writer.flush();
        writer.close();
    }

    public static void main(String[] args) throws Exception {
        //generateJunkIndex();
    }
}
