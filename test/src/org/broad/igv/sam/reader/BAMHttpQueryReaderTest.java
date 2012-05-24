/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.sam.reader;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.sam.Alignment;
import org.broad.igv.util.ResourceLocator;
import org.junit.*;

import java.io.File;
import java.io.IOException;

import static org.junit.Assert.*;

/**
 * @author jrobinso
 */
public class BAMHttpQueryReaderTest extends AbstractHeadlessTest {

    //private final String BAM_URL_STRING = "http://www.broadinstitute.org/igvdata/test/index_test.bam";
    private final String BAM_URL_STRING = "http://www.broadinstitute.org/igvdata/1KG/freeze5_merged/low_coverage_CEU.Y.bam";

    BAMHttpReader reader;

    @BeforeClass
    public static void setUpClass() throws Exception {
        AbstractHeadlessTest.setUpClass();
        SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.SILENT);
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
        SAMFileHeader header = reader.getHeader();
        assertEquals(114, header.getSequenceDictionary().size());
        assertEquals("1.0", header.getVersion());
    }

    @Test
    public void testIterator() {
        CloseableIterator<Alignment> iter = reader.iterator();
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
        CloseableIterator<Alignment> iter = reader.query(chr, start, end, false);
        int counted = 0;
        while (iter.hasNext()) {
            Alignment a = iter.next();
            counted++;
            assertNotNull(a);
        }
        iter.close();

        assertEquals(expected_count, counted);
    }

    @Test
    public void testDeleteIndex() throws Exception {
        File indexFile = reader.indexFile;

        assertTrue(indexFile.exists());

        CloseableIterator<Alignment> iter = reader.query("Y", 10000000, 10004000, false);
        int max = 100;
        int counted = 0;
        while (iter.hasNext()) {
            Alignment a = iter.next();
            counted++;
            assertNotNull(a);
            if (counted > max) {
                break;
            }
        }

        //Not closing the iterator on purpose

        assertTrue(indexFile.exists());

        assertTrue(indexFile.canWrite());
    }

}