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

package org.broad.igv.bbfile;

import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.*;

/**
 * @author jrobinso
 * @date Jun 19, 2011
 */
public class BBFileReaderTest {


    @Test
    public void testBigBed() throws IOException {

        String path = TestUtils.DATA_DIR + "bb/chr21.refseq.bb";
        BBFileReader bbReader = new BBFileReader(path);

        BBFileHeader bbFileHdr = bbReader.getBBFileHeader();
        assertTrue(bbFileHdr.isBigBed());

        String chr = "chr21";
        int start = 26490012;
        int end = 42182827;


        for (BBZoomLevelHeader header : bbReader.getZoomLevels().getZoomLevelHeaders()) {
            assertNotNull(header);

            ZoomLevelIterator zlIter = bbReader.getZoomLevelIterator(header.getZoomLevel(), chr, start, chr, end, false);

            while (zlIter.hasNext()) {
                ZoomDataRecord rec = zlIter.next();
                int n = rec.getBasesCovered();
                if (n > 0) {
                    assertEquals(chr, rec.getChromName());
                    assertTrue(rec.getChromEnd() >= start && rec.getChromStart() <= end);
                }

            }
        }


    }


}
