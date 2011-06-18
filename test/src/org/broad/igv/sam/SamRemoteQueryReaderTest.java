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
package org.broad.igv.sam;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.sam.reader.BAMRemoteQueryReader;
import org.broad.igv.util.ResourceLocator;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;


/**
 * @author jrobinso
 */
public class SamRemoteQueryReaderTest {

    public SamRemoteQueryReaderTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    /**
     * Test of query method, of class SamRemoteQueryReader.
     */
    @Test
    public void testQuery() {
        String path = "/broad/1KG/DCC_merged/freeze4/NA12878.ceu.daughter.bam";
        String serverURL = "http://node253:7080/webservices/igv";

        String chr = "1";
        int start = 713700;
        int end = 714100;

        BAMRemoteQueryReader reader = new BAMRemoteQueryReader(new ResourceLocator(serverURL, path));

        SAMFileHeader header = reader.getHeader();
        System.out.println(header.getSequence("1"));

        CloseableIterator<Alignment> iter = reader.query(chr, start, end, true);

        //SAMRecord firstRecord = iter.next();
        //assertEquals(expectedFirstRecord, firstRecord.format());
        while (iter.hasNext()) {
            Alignment record = iter.next();
            System.out.println(record.getReadName());
        }

        iter.close();

    }
}