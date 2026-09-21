
package org.igv.alignment;

import htsjdk.samtools.util.CloseableIterator;
import org.igv.Globals;
import org.igv.alignment.reader.BAMReader;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

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

    @Test
    public void testCSI() throws Exception {

        String bamfile = TestUtils.DATA_DIR + "bam/BAMFileIndexTest/index_test.bam";
        String baifile = TestUtils.DATA_DIR + "bam/BAMFileIndexTest/index_test.bam.bai";
        String csifile = TestUtils.DATA_DIR + "bam/BAMFileIndexTest/index_test.bam.csi";

        String chr = "chr1";
        int end = 6000000;
        int start = 1000000;

        ResourceLocator baiLocator = new ResourceLocator(bamfile);
        baiLocator.setIndexPath(baifile);
        BAMReader baireader = new BAMReader(baiLocator, true);

        ResourceLocator csiLocator = new ResourceLocator(bamfile);
        csiLocator.setIndexPath(csifile);
        BAMReader csiReader = new BAMReader(csiLocator, true);



        CloseableIterator<SAMAlignment> baiiter = baireader.query(chr, start, end, true);
        CloseableIterator<SAMAlignment> csiiter = csiReader.query(chr, start, end, true);

        int count = 0;
        while (baiiter.hasNext()) {
            Alignment bamrecord = baiiter.next();
            Alignment samrecord = csiiter.next();
            assertTrue(bamrecord.getStart() >= start);
            assertTrue(bamrecord.getEnd() <= end);
            assertEquals(bamrecord.getReadName(), samrecord.getReadName());
            assertEquals(bamrecord.getSample(), samrecord.getSample());
            count++;
        }
        assertTrue("Unexpected data count: " + count, count == 20);

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