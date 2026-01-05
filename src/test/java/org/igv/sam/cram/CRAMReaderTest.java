
package org.igv.sam.cram;

import htsjdk.samtools.util.CloseableIterator;
import org.igv.Globals;
import org.igv.feature.genome.GenomeManager;
import org.igv.sam.SAMAlignment;
import org.igv.sam.reader.BAMReader;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * @author jrobinso
 */
public class CRAMReaderTest {

    public CRAMReaderTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        Globals.setHeadless(true);
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }


    @Test
    public void testIterateLocalCraiCram() throws Exception {

        String cramFile = TestUtils.DATA_DIR + "cram/cram_with_crai_index.cram";
        GenomeManager.getInstance().loadGenome(TestUtils.DATA_DIR + "cram/hg19mini.fasta");

        BAMReader reader = new BAMReader(new ResourceLocator(cramFile), true);

        List<String> seqNames = reader.getSequenceNames();
        assertEquals(4, seqNames.size());


        CloseableIterator<SAMAlignment> iter = reader.iterator();
        int counter = count(iter);
        assertEquals(11, counter);

    }

    @Test
    public void testQueryLocalCraiCram() throws Exception {

        String cramFile = TestUtils.DATA_DIR + "cram/cram_with_crai_index.cram";
        GenomeManager.getInstance().loadGenome(TestUtils.DATA_DIR + "cram/hg19mini.fasta");

        BAMReader reader = new BAMReader(new ResourceLocator(cramFile), true);

        CloseableIterator<SAMAlignment> iter = reader.query("2", 500, 600, false);
        int counter = count(iter);
        assertEquals(2, counter);
    }


    @Test
    public void testQueryLocalBaiCram() throws Exception {

        String cramFile = TestUtils.DATA_DIR + "cram/cram_with_bai_index.cram";
        GenomeManager.getInstance().loadGenome(TestUtils.DATA_DIR + "cram/hg19mini.fasta");

        BAMReader reader = new BAMReader(new ResourceLocator(cramFile), true);

        CloseableIterator<SAMAlignment> iter = reader.query("2", 500, 600, false);
        int counter = count(iter);
        assertEquals(2, counter);
    }

    @Test
    public void testRemoteBaiCram() throws Exception {

        String cramFile = TestUtils.DATA_DIR + "cram/cram_with_bai_index.cram";
        GenomeManager.getInstance().loadGenome(TestUtils.DATA_DIR + "cram/hg19mini.fasta");

        BAMReader reader = new BAMReader(new ResourceLocator(cramFile), true);

        CloseableIterator<SAMAlignment> iter = reader.query("2", 500, 600, false);
        int counter = count(iter);
        assertEquals(2, counter);
    }


    public int count(CloseableIterator<SAMAlignment> iter) {
        int counter = 0;
        while (iter.hasNext()) {
            iter.next();
            counter++;
        }
        return counter;
    }

}