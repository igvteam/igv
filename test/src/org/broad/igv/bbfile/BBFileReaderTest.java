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

/**
 * @author jrobinso
 * @date Jun 19, 2011
 */
public class BBFileReaderTest {


    public static BBFileReader openReader(String path) throws IOException {

        // time the B+/R+ chromosome an zoom level tree construction
        long time = System.currentTimeMillis(), time_prev = time;

        BBFileReader bbReader = new BBFileReader(path);

        // get the time mark
        time = System.currentTimeMillis();
        System.out.println("B+/R+/Zoom tree build = " + (time - time_prev) + " ms");

        // check file type
        BBFileHeader bbFileHdr = bbReader.getBBFileHeader();
        if (bbReader.isBigBedFile())
            System.out.println("BBFileReader header indicates a BigBed file type\n");
        else if (bbReader.isBigWigFile())
            System.out.println("BBFileReader header indicates a BigWig file type\n");
        else
            System.out.println("BBFileReader was not able to identify file type\n");

        return bbReader;
    }

    @Test
    public void testBigBed() throws IOException {

        String path = TestUtils.DATA_DIR + "/bb/chr21.refseq.bb";
        BBFileReader bbReader = openReader(path);

        String chr = "chr21";
        int start = 26490012;
        int end = 42182827;


        for (BBZoomLevelHeader header : bbReader.getZoomLevels().getZoomLevelHeaders()) {

            header.print();
        }


        ///end = 100000;
        ZoomLevelIterator zlIter = bbReader.getZoomLevelIterator(3, chr, start, chr, end, false);

        while (zlIter.hasNext()) {
            ZoomDataRecord rec = zlIter.next();
            int n = rec.getBasesCovered();
            if (n > 0) {
                System.out.println("Base covered = " + n);
                double mean = rec.getSumData() / n;
                System.out.println(rec.getChromName() + "\t" + rec.getChromStart() + "\t" + rec.getChromEnd() +
                        "\t" + mean);
            }

        }


    }

}
