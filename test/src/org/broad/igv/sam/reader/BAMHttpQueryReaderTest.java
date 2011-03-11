/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.CachingQueryReader;
import org.junit.AfterClass;
import static org.junit.Assert.assertEquals;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.net.URL;

/**
 * @author jrobinso
 */
public class BAMHttpQueryReaderTest {

    //private final String BAM_URL_STRING = "http://www.broadinstitute.org/igvdata/test/index_test.bam";
    private final String BAM_URL_STRING = "http://www.broadinstitute.org/igvdata/1KG/DCC_merged/freeze5/NA12891.pilot2.SLX.bam";

    AlignmentQueryReader reader;

    public BAMHttpQueryReaderTest() {


    }

    @Before
    public void setUpClass() throws Exception {
        URL url = new URL(BAM_URL_STRING);
        CachingQueryReader reader = new CachingQueryReader(new BAMHttpQueryReader(url, true));
        System.out.println("Index loaded");
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Test
    public void testClose() throws Exception {
    }

    @Test
    public void testGetHeader() throws IOException {
        SAMFileHeader header = reader.getHeader();
        assertEquals(45, header.getSequenceDictionary().size());

    }

    @Test
    public void testIterator() {
        CloseableIterator<Alignment> iter = reader.iterator();
        while (iter.hasNext()) {
            Alignment a = iter.next();
            System.out.println(a.getReadName());
        }

    }

    @Test
    public void testQuery() throws IOException {
        SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        CloseableIterator<Alignment> iter = reader.query("1", 100000000, 100004000, false);
        while (iter.hasNext()) {
            Alignment a = iter.next();
            System.out.println(a.getReadName());

        }
    }

}