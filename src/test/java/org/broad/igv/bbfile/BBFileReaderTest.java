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
