/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.igv.sam.reader;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.CloseableIterator;
import org.igv.AbstractHeadlessTest;
import org.igv.sam.Alignment;
import org.igv.sam.SAMAlignment;
import org.igv.util.ResourceLocator;
import org.junit.*;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import static org.junit.Assert.*;

/**
 * @author jrobinso
 */

public class BAMHttpReaderTest extends AbstractHeadlessTest {

    private static final String BAM_URL_STRING = "https://1000genomes.s3.amazonaws.com/phase3/data/HG01879/exome_alignment/HG01879.mapped.ILLUMINA.bwa.ACB.exome.20120522.bam";

    BAMReader reader;

    @BeforeClass
    public static void setUpClass() throws Exception {
        AbstractHeadlessTest.setUpClass();
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() throws Exception {
        reader = new BAMReader(new ResourceLocator(BAM_URL_STRING), true);
    }

    @After
    public void tearDown() throws Exception {
        if (reader != null) reader.close();
        reader = null;
    }

    @Test
    public void testGetHeader() throws IOException {
        SAMFileHeader header = reader.getFileHeader();
        assertEquals(86, header.getSequenceDictionary().size());
        assertEquals("1.0", header.getVersion());
    }

    @Test
    public void testIterator() throws IOException {
        CloseableIterator<SAMAlignment> iter = reader.iterator();
        //This will iterate over the entire file, so we break and exit after a few iterations
        int minnum = 10;
        int actnum = 0;
        while (iter.hasNext()) {
            Alignment a = iter.next();
            assertNotNull(a);
            actnum++;

            if (actnum > minnum) {
                break;
            }
        }
        iter.close();
        assertTrue(actnum > minnum);

    }

    @Test
    public void testQuery() throws Exception {
        int expected_count = 4;
        String chr = "Y";
        int start = 10000000 - 1;
        int end = 10004000;
        
        CloseableIterator<SAMAlignment> iter = reader.query(chr, start, end, false);
        int counted = 0;
        while (iter.hasNext()) {
            Alignment a = iter.next();
            counted++;
            assertNotNull(a);
        }
        iter.close();

        assertEquals(expected_count, counted);
    }


}
