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
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.PicardAlignment;
import org.broad.igv.util.ResourceLocator;
import org.junit.*;

import java.io.File;
import java.io.IOException;

import static org.junit.Assert.*;

/**
 * @author jrobinso
 */
@Ignore
public class BAMHttpQueryReaderTest extends AbstractHeadlessTest {

    //private final String BAM_URL_STRING = "http://www.broadinstitute.org/igvdata/test/index_test.bam";
    private final String BAM_URL_STRING = "http://www.broadinstitute.org/igvdata/1KG/freeze5_merged/low_coverage_CEU.Y.bam";

    BAMHttpReader reader;

    @BeforeClass
    public static void setUpClass() throws Exception {
        AbstractHeadlessTest.setUpClass();
        SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() throws Exception {
        reader = new BAMHttpReader(new ResourceLocator(BAM_URL_STRING), true);
    }

    @After
    public void tearDown() throws Exception {
        reader.close();
        reader = null;
    }

    @Test
    public void testGetHeader() throws IOException {
        SAMFileHeader header = reader.getFileHeader();
        assertEquals(114, header.getSequenceDictionary().size());
        assertEquals("1.0", header.getVersion());
    }

    @Test
    public void testIterator() {
        CloseableIterator<PicardAlignment> iter = reader.iterator();
        //This takes a long time. We just look for a minimum number
        int minnum = 1000000;
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

        checkNumber("Y", 10000000 - 1, 10004000, 6890);
        checkNumber("Y", 100000000 - 1, 100040000, 0);
        checkNumber("1", 1 - 1, 100000000, 0);

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