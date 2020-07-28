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

package org.broad.igv.feature.genome;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.DirectoryManager;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.commandbar.GenomeListManager;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.TestUtils;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.net.URL;

import static org.junit.Assert.*;

/**
 * @author jrobinso
 */
public class GenomeManagerTest extends AbstractHeadlessTest {

    static GenomeManager genomeManager;

    public GenomeManagerTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        GenomeManager.getInstance().clearGenomeCache();
        AbstractHeadlessTest.setUpClass();
        genomeManager = GenomeManager.getInstance();
    }

    @Test
    public void testGenerateGenomeList() throws Exception {
        File inDir = new File(TestUtils.DATA_DIR, "genomes");
        String outPath = TestUtils.TMP_OUTPUT_DIR + "/genomelist.txt";

        String rootPath = "http://igvdata.broadinstitute.org/genomes";

        genomeManager.generateGenomeList(inDir, rootPath, outPath);

        BufferedReader reader = new BufferedReader(new FileReader(outPath));
        int count = 0;
        String line;
        while ((line = reader.readLine()) != null) {
            assertTrue(line.contains(rootPath + "/"));
            count++;
        }
        assertEquals(6, count);
    }

    private String genomeZipFile = TestUtils.TMP_OUTPUT_DIR + "tmp.genome";
    private String fastaFileRelPath = TestUtils.DATA_DIR + "fasta/ecoli_out.padded.fasta";
    private String genomeDisplayName = "Unit test genome";
    private String genomeId = "gmt_001";

    private void createDotGenomeForTest(String fastaFileName) throws IOException{
        GenomeListManager.getInstance().getUserDefinedGenomeMap();
        GenomeListItem genomeListItem = GenomeManager.getInstance().defineGenome(
                new File(genomeZipFile), null, null,
                fastaFileName, null, genomeDisplayName,
                genomeId, null);


    }

    /**
     * Use a relative path for fasta file, test that we can load it
     * @throws Exception
     */
    @Test
    public void testLoadGenomeFastaRelative() throws Exception{
        createDotGenomeForTest(fastaFileRelPath);
        Genome relGenome = GenomeManager.getInstance().loadGenome(genomeZipFile, null);

        checkGenome(relGenome);
    }

    /**
     * Use an absolute path for fasta file, test that we can load it
     * @throws Exception
     */
    @Test
    public void testLoadGenomeFastaAbsolute() throws Exception{
        File fastaFile = new File(fastaFileRelPath);
        String fastaAbsPath = fastaFile.getAbsolutePath();

        createDotGenomeForTest(fastaAbsPath);
        Genome absGenome = GenomeManager.getInstance().loadGenome(genomeZipFile, null);

        checkGenome(absGenome);
    }

    private void checkGenome(Genome genome) {
        String chr = genome.getAllChromosomeNames().get(0);
        int end = 10;

        byte[] seq = genome.getSequence(chr, 0, end);
        assertNotNull(seq);
        assertEquals(end, seq.length);
    }

    @Test
    public void testLoadFastaOrdering() throws Exception{
        String fastaPath = TestUtils.DATA_DIR + "fasta/out_order.fa";
        TestUtils.createIndex(fastaPath);

        Genome genome = GenomeManager.getInstance().loadGenome(fastaPath, null);
        String[] chromos = {"chr5", "chr1"};

        assertArrayEquals(chromos, genome.getAllChromosomeNames().toArray());
    }

    /**
     * Test defining a genome from a chrom.sizes file
     *
     * @throws Exception
     */
    @Test
    public void testLoadChromSizes() throws Exception {
        String testFile = TestUtils.DATA_DIR + "genomes/hg19.chrom.sizes";
        Genome genome = GenomeManager.getInstance().loadGenome(testFile, null);

        assertEquals(24, genome.getAllChromosomeNames().size());
        assertEquals(3095677412l, genome.getTotalLength());
    }

}
