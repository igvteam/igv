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
import org.broad.igv.PreferenceManager;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.RunnableResult;
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
        GenomeManager.getInstance().getUserDefinedGenomeArchiveList();
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

    /**
     * Test for simple downloading of .genome file,
     * which packages up all the resources it can
     * @throws Exception
     */
    @Test
    public void testDownloadDotGenome() throws Exception{

        String genId = "NC_001802";
        String genomePath = "http://igv.broadinstitute.org/genomes/" + genId + ".genome";
        String outDirPath = TestUtils.TMP_OUTPUT_DIR;
        File outDirFile = new File(outDirPath);
        File outGenomeFile = new File(outDirPath, genId + ".genome");
        RunnableResult result = GenomeManager.getInstance().downloadWholeGenome(genomePath, outDirFile, null);
        assertTrue("Download of genome failed", result.isSuccess());

        assertTrue(outGenomeFile.exists());
        File fastaFile = new File(TestUtils.TMP_OUTPUT_DIR, genId + ".fna");
        assertTrue("fasta file not found: " + fastaFile.getAbsolutePath(), fastaFile.exists());

        //I don't know exactly why this is necessary but it seems like there is some weird memory-buffering issue
        //and if we try to read from the outPath again it uses the old (remote) fasta file name.
        //This is just a unit test thing, copying the downloaded .genome file should have the same info
        //and it makes the tests pass
        File tmpOut = new File(TestUtils.TMP_OUTPUT_DIR + "t2.genome");
        FileUtils.copyFile(outGenomeFile, tmpOut);
        GenomeDescriptor descriptor = GenomeManager.parseGenomeArchiveFile(tmpOut);
        assertEquals(fastaFile.getAbsolutePath(), descriptor.getSequenceLocation());
        assertTrue(descriptor.hasCustomSequenceLocation());

        String remSequencePath = "http://igvdata.broadinstitute.org/genomes/seq/Human_immunodeficiency_virus_1_uid15476/NC_001802.fna";
        FastaIndexedSequence remSequence = new FastaIndexedSequence(remSequencePath);

        TestUtils.createIndex(fastaFile.getAbsolutePath());
        FastaIndexedSequence localSequence = new FastaIndexedSequence(fastaFile.getAbsolutePath());

        assertEquals(remSequence.getChromosomeNames(), localSequence.getChromosomeNames());
    }

    @Test
    public void testRewriteSequenceLocation() throws Exception{
        String origGenomePath = TestUtils.DATA_DIR + "genomes/hg18.unittest.genome";
        File newGenomeFile = File.createTempFile("hg18.unittest.testrewrite", ".genome");
        FileUtils.copyFile(new File(origGenomePath), newGenomeFile);
        newGenomeFile.deleteOnExit();

        String newSeqLocation = (new File(TestUtils.TMP_OUTPUT_DIR + "/myseq.fa")).getAbsolutePath();
        boolean success = GenomeManager.rewriteSequenceLocation(newGenomeFile, newSeqLocation);

        assertTrue("Rewrite of .genome failed", success);
        GenomeDescriptor descriptor = GenomeManager.parseGenomeArchiveFile(newGenomeFile);

        assertEquals(newSeqLocation, descriptor.getSequenceLocation());
        assertTrue(descriptor.hasCustomSequenceLocation());
    }

    @Test
    public void testLoadNonCustomGenome() throws Exception{
        String genomePath = TestUtils.DATA_DIR + "genomes/hg18.unittest.genome";
        GenomeDescriptor descriptor = GenomeManager.parseGenomeArchiveFile(new File(genomePath));
        assertFalse(descriptor.hasCustomSequenceLocation());
    }

    /**
     * Test that when we update a cached genome with custom sequence, the sequence location is preserved
     * We have 2 .genome files in the test directory which are very similar.
     * We copy one to the genome cache directory and rewrite the sequence
     * Then we download the other (needs to be remote)
     * We then check that the rewritten sequence was preserved, and differing properties were updated
     * @throws Exception
     */

    // Test disabled -- unreliable on the server
    //@Test
    public void testRefreshCacheLocalSequence() throws Exception{
        String cachedSrcPath = TestUtils.DATA_DIR + "genomes/local.unittest_cached.genome";
        File cachedSrcFile = new File(cachedSrcPath);
        File cachedFile = new File(DirectoryManager.getGenomeCacheDirectory(), cachedSrcFile.getName());
        FileUtils.copyFile(cachedSrcFile, cachedFile);

        //The local variable names "updatedDescriptor" and "newDescriptor" are somewhat
        //confusing, but I don't know how to do better
        File updatedFile = new File(TestUtils.DATA_DIR + "genomes/local.unittest_updated.genome");
        GenomeDescriptor updatedDescriptor = GenomeManager.parseGenomeArchiveFile(updatedFile);

        GenomeDescriptor cachedDescriptor = GenomeManager.parseGenomeArchiveFile(cachedFile);
        String targetSeqLocation = cachedDescriptor.getSequenceLocation();
        assert (cachedDescriptor.hasCustomSequenceLocation());

        //TODO We could point to github instead of the Broad, only real benefit is if we ever change
        //this file it's one less manual update step
        //String parent = "https://github.com/broadinstitute/IGV/tree/master/test/data/genomes";
        String parent = "http://www.broadinstitute.org/igvdata/test";
        URL remURL = new URL(parent + "/" + updatedFile.getName());

        PreferenceManager.getInstance().put(PreferenceManager.AUTO_UPDATE_GENOMES, true);
        GenomeManager.getInstance().refreshCache(cachedFile, remURL);

        GenomeDescriptor newDescriptor = GenomeManager.parseGenomeArchiveFile(cachedFile);
        assertEquals(targetSeqLocation, newDescriptor.getSequenceLocation());
        assertTrue(newDescriptor.hasCustomSequenceLocation());

        assertEquals(updatedDescriptor.getName(), newDescriptor.getName());
        assertEquals(updatedDescriptor.getGeneFileName(), newDescriptor.getGeneFileName());

        assertNotSame(newDescriptor.getName(), cachedDescriptor.getName());
    }



}