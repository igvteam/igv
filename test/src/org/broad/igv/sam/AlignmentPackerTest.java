/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.sam;

import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.List;
import java.util.Map;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * @author jrobinso
 */
public class AlignmentPackerTest {

    public AlignmentPackerTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        TestUtils.setUpHeadless();
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    /**
     * Test of packAlignments method, of class AlignmentPacker.
     */
    @Test
    public void testPackAlignments() throws Exception {

        //Retrieve data to pack
        String path = "http://www.broadinstitute.org/igvdata/1KG/pilot2Bams/NA12878.SLX.bam";
        String chr = "1";
        int start = 557000;
        int end = 558000;
        boolean contained = false;

        ResourceLocator rl = new ResourceLocator(path);
        AlignmentReader samReader = AlignmentReaderFactory.getReader(rl);
        CloseableIterator<Alignment> iter = samReader.query(chr, start, end, contained);
        ///////////////////////////

        boolean showDuplicates = false;
        int qualityThreshold = 0;
        int maxLevels = 1000;

        Map<String, List<AlignmentInterval.Row>> result = (new AlignmentPacker()).packAlignments(iter, end, false, new AlignmentTrack.RenderOptions(), 10000);
        assertEquals(1, result.size());
        for (List<AlignmentInterval.Row> alignmentrows : result.values()) {
            for (AlignmentInterval.Row alignmentrow : alignmentrows) {
                List<Alignment> alignments = alignmentrow.alignments;
                for (int ii = 1; ii < alignments.size(); ii++) {
                    assertTrue(alignments.get(ii).getAlignmentStart() > alignments.get(ii - 1).getAlignmentStart());
                    assertTrue(alignments.get(ii).getAlignmentStart() - alignments.get(ii - 1).getAlignmentEnd() >= AlignmentPacker.MIN_ALIGNMENT_SPACING);
                }
            }
        }


    }


}