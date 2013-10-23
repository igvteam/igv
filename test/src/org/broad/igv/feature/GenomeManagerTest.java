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

package org.broad.igv.feature;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.genome.*;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.RunnableResult;
import org.broad.igv.util.TestUtils;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

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
        String outPath = TestUtils.DATA_DIR + "out/genomelist.txt";

        String rootPath = "http://igvdata.broadinstitute.org/genomes";

        genomeManager.generateGenomeList(inDir, rootPath, outPath);

        BufferedReader reader = new BufferedReader(new FileReader(outPath));
        int count = 0;
        String line;
        while ((line = reader.readLine()) != null) {
            assertTrue(line.contains(rootPath + "/"));
            count++;
        }
        assertEquals(5, count);
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
        RunnableResult result = GenomeManager.getInstance().downloadWholeGenome(genomePath, outDirFile);
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
        GenomeDescriptor descriptor = GenomeManager.getInstance().parseGenomeArchiveFile(tmpOut);
        assertEquals(fastaFile.getAbsolutePath(), descriptor.getSequenceLocation());

        String remSequencePath = "http://igvdata.broadinstitute.org/genomes/seq/Human_immunodeficiency_virus_1_uid15476/NC_001802.fna";
        FastaIndexedSequence remSequence = new FastaIndexedSequence(remSequencePath);

        TestUtils.createIndex(fastaFile.getAbsolutePath());
        FastaIndexedSequence localSequence = new FastaIndexedSequence(fastaFile.getAbsolutePath());

        assertEquals(remSequence.getChromosomeNames(), localSequence.getChromosomeNames());
    }



}