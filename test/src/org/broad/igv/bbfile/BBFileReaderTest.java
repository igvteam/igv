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

import org.junit.BeforeClass;
import org.junit.Test;

import java.io.IOException;
import java.util.Iterator;

/**
 * @author jrobinso
 * @date Jun 19, 2011
 */
public class BBFileReaderTest {


    public static BBFileReader openReader(String path) throws IOException {

        // time the B+/R+ chromosome an zoom level tree construction
        long time = System.currentTimeMillis(), time_prev = time;

        BBFileReader bbReader = new BBFileReader("/Users/jrobinso/projects/bigwig/test/data/ce6_uniqueome.coverage.base-space.25.1.bw");

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

        String path = "/Users/jrobinso/projects/bigwig/test/data/chr21.bb";
        BBFileReader bbReader = openReader(path);

        //chr21:26,490,012-42,182,827

        String chr = "chr21";
        int start = 26490012;
        int end = 42182827;


        for (BBZoomLevelHeader header : bbReader.getZoomLevels().getZoomLevelHeaders()) {

            header.print();
        }


        ///end = 100000;
        ZoomLevelIterator zlIter = bbReader.getZoomLevelIterator(7, chr, start, chr, end, false);

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

    //@Test
    public void testBigWig() throws IOException {

        String path = "/Users/jrobinso/projects/bigwig/test/data/ce6_uniqueome.coverage.base-space.25.1.bw";
        BBFileReader bbReader = openReader(path);

        String chr = "chrI";
        int start = 7517811;
        int end = 7521491;


        ///end = 100000;
        ZoomLevelIterator zlIter = bbReader.getZoomLevelIterator(1, chr, start, chr, end, false);

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


        Iterator<WigItem> iter = bbReader.getBigWigIterator(chr, start, chr, end, false);

        double sum = 0;
        int nb = 0;
        int basesCovered = 0;
        while (iter.hasNext()) {
            WigItem wi = iter.next();
            //System.out.println(wi.getChromosome() + "\t" + wi.getStartBase() + "\t" + wi.getEndBase()
            //        + "\t" + wi.getWigValue());

            double span = wi.getEndBase() - wi.getStartBase();
            sum += span * wi.getWigValue();
            nb++;
            basesCovered += span;
        }
        System.out.println("mean = " + (sum / basesCovered));
        System.out.println("nb = " + nb);
        System.out.println("basesCovere = " + basesCovered);


        System.out.println("Zoom levels");
        //BBZoomLevels levels = bbReader.getZoomLevels();

        for (BBZoomLevelHeader header : bbReader.getZoomLevels().getZoomLevelHeaders()) {
            header.print();
        }

    }
}
