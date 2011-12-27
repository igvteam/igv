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

import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.sam.reader.BAMRemoteQueryReader;
import org.broad.igv.util.ResourceLocator;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
public class AlignmentPackerTest {

    public AlignmentPackerTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    /**
     * Test of packAlignments method, of class AlignmentPacker.
     */
    @Test
    public void testPackAlignments() {
        /*
        String path = "/Users/jrobinso/IGV/Alignments/303KY.8.paired.bam";
        String serverURL = "http://localhost:8080/webservices/igv";
        String chr = "chr7";
        int start = 2244149;
        int end = 2250144;
         * */

        String path = "http://www.broadinstitute.org/igvdata/1KG/pilot2Bams/NA12878.SLX.bam";
        String chr = "1";
        int start = 557000;  //98751004; //
        int end = 558000; //98751046; //
        boolean contained = false;

        ResourceLocator rl = new ResourceLocator(path);

        BAMRemoteQueryReader bamReader = new BAMRemoteQueryReader(rl);
        CloseableIterator<Alignment> iter = bamReader.query(chr, start, end, contained);
        boolean showDuplicates = false;
        int qualityThreshold = 0;
        int maxLevels = 1000;


        Map<String, List<AlignmentInterval.Row>> result = (new AlignmentPacker()).packAlignments(iter, end, false, null, 10000);

    }


}