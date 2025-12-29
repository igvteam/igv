/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.igv.feature.genome;

import org.igv.AbstractHeadlessTest;
import org.igv.util.TestUtils;
import org.junit.BeforeClass;
import org.junit.Test;

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
        AbstractHeadlessTest.setUpClass();
        genomeManager = GenomeManager.getInstance();
    }


    @Test
    public void testLoadFastaOrdering() throws Exception{
        String fastaPath = TestUtils.DATA_DIR + "fasta/out_order.fa";
        TestUtils.createIndex(fastaPath);

        Genome genome = GenomeManager.getInstance().loadGenome(fastaPath);
        String[] chromos = {"chr5", "chr1"};

        assertArrayEquals(chromos, genome.getChromosomeNames().toArray());
    }

    /**
     * Test defining a genome from a chrom.sizes file
     *
     * @throws Exception
     */
    @Test
    public void testLoadChromSizes() throws Exception {
        String testFile = TestUtils.DATA_DIR + "genomes/hg19.chrom.sizes";
        Genome genome = GenomeManager.getInstance().loadGenome(testFile);

        assertEquals(37, genome.getChromosomeNames().size());
    }

}
