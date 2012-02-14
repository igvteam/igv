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

package org.broad.igv.bbfile;

import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

/**
 * @author jrobinso
 * @date Jun 19, 2011
 */
public class BBFileReaderTest {


    @Test
    public void testBigBed() throws IOException {

        String path = TestUtils.DATA_DIR + "/bb/chr21.refseq.bb";
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
