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

import org.broad.igv.bigwig.BBDataSource;
import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.*;

/**
 * @author jrobinso
 * @date Jun 19, 2011
 */
public class BWSourceTest {


    @Test
    public void testBigGenePred() throws IOException {

        String path = "https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/bbi/GCA_009914755.4_T2T-CHM13v2.0.catLiftOffGenesV1/catLiftOffGenesV1.bb";

        BBFileReader bbReader = new BBFileReader(path);
        BBDataSource bbSource = new BBDataSource(bbReader, null);

        BBFileHeader bbFileHdr = bbReader.getBBFileHeader();
        assertTrue(bbFileHdr.isBigBed());

        String autoSql = bbReader.getAutoSql();
System.out.println(autoSql);
        String chr = "chr21";
        int start = 26490012;
        int end = 42182827;


        BigBedIterator iter = bbReader.getBigBedIterator(chr, start, chr, end, false);
        int count = 0;
        while (iter.hasNext()) {
            BedFeature f = iter.next();
            assertEquals(chr, f.getChromosome());
            assertTrue(f.getStartBase() <= end && f.getEndBase() >= start);
            count++;
        }
        assertEquals("Feature count", 225, count);

    }


}
