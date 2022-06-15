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

import htsjdk.tribble.Feature;
import org.junit.Test;

import java.io.IOException;
import java.util.Iterator;

import static org.junit.Assert.*;

/**
 * @author jrobinso
 * @date Jun 19, 2011
 */
public class BBDataSourceTest {


    @Test
    public void testBigGenePred() throws IOException {

        String path = "https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/bbi/GCA_009914755.4_T2T-CHM13v2.0.catLiftOffGenesV1/catLiftOffGenesV1.bb";

        BBFileReader bbReader = new BBFileReader(path);
        BBDataSource bbSource = new BBDataSource(bbReader, null);

        String chr = "chr21";
        int start = 26490012;
        int end = 42182827;

        Iterator<Feature> iter = bbSource.getFeatures(chr, start, end);

        int count = 0;
        while (iter.hasNext()) {
            Feature f = iter.next();
            assertEquals(chr, f.getChr());
            assertTrue(f.getStart() <= end && f.getEnd() >= start);
            count++;
        }
        assertTrue(count > 0);

    }

    @Test
    public void testBigNarrowpeak() throws IOException {

        String path = "https://www.encodeproject.org/files/ENCFF154CVA/@@download/ENCFF154CVA.bigBed";

        BBFileReader bbReader = new BBFileReader(path);
        BBDataSource bbSource = new BBDataSource(bbReader, null);

        String chr = "chr21";
        int start = 26490012;
        int end = 42182827;

        Iterator<Feature> iter = bbSource.getFeatures(chr, start, end);

        int count = 0;
        while (iter.hasNext()) {
            Feature f = iter.next();
            assertEquals(chr, f.getChr());
            assertTrue(f.getStart() <= end && f.getEnd() >= start);
            count++;
        }
        assertTrue(count > 0);

    }

    // NOTE:  bigInteract is not yet supported in IGV desktop, bigInteract is treated as a plain bed file.
    @Test
    public void testBigInteract() throws IOException {

        String path = "https://s3.amazonaws.com/igv.org.demo/interactExample3.inter.bb";

        BBFileReader bbReader = new BBFileReader(path);
        BBDataSource bbSource = new BBDataSource(bbReader, null);

        String chr = "chr3";
        int start = 63081504;
        int end = 64501215;

        Iterator<Feature> iter = bbSource.getFeatures(chr, start, end);

        int count = 0;
        while (iter.hasNext()) {
            Feature f = iter.next();
            assertEquals(chr, f.getChr());
            assertTrue(f.getStart() <= end && f.getEnd() >= start);
            count++;
        }
        assertTrue(count > 0);
    }

    // NOTE:  bigInteract is not yet supported in IGV desktop, bigInteract is treated as a plain bed file.
    @Test
    public void testBigRmsk() throws IOException {

        String path = "https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/bbi/GCA_009914755.4_T2T-CHM13v2.0.t2tRepeatMasker/chm13v2.0_rmsk.bb";

        BBFileReader bbReader = new BBFileReader(path);
        BBDataSource bbSource = new BBDataSource(bbReader, null);

        String chr = "chr9";
        int start = 145481787;
        int end = 145481785;

        Iterator<Feature> iter = bbSource.getFeatures(chr, start, end);

        int count = 0;
        while (iter.hasNext()) {
            Feature f = iter.next();
            assertEquals(chr, f.getChr());
            assertTrue(f.getStart() <= end && f.getEnd() >= start);
            count++;
        }
        assertTrue(count > 0);
    }

}
