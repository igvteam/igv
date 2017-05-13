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

package org.broad.igv.sam.reader;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.CloseableIterator;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.PicardAlignment;
import org.broad.igv.util.ResourceLocator;
import org.junit.*;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import static org.junit.Assert.*;

/**
 * @author jrobinso
 */

public class BAMHttpReaderTest extends AbstractHeadlessTest {

    private static final String BAM_URL_STRING = "http://1000genomes.s3.amazonaws.com/phase3/data/HG01879/exome_alignment/HG01879.mapped.ILLUMINA.bwa.ACB.exome.20120522.bam";

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
        reader.close();
        reader = null;
    }

    @Test
    public void testGetHeader() throws IOException {
        SAMFileHeader header = reader.getFileHeader();
        assertEquals(86, header.getSequenceDictionary().size());
        assertEquals("1.0", header.getVersion());
    }

    @Test
    public void testIterator() {
        CloseableIterator<PicardAlignment> iter = reader.iterator();
        //This takes a long time. We just look for a minimum number
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
        checkNumber("Y", 10000000 - 1, 10004000, 4);
    }

    private void checkNumber(String chr, int start, int end, int expected_count) throws IOException {
        CloseableIterator<PicardAlignment> iter = reader.query(chr, start, end, false);
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
