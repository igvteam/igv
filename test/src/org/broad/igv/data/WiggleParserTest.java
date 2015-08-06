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

package org.broad.igv.data;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

/**
 * @author jacob
 * @date 2013-Aug-12
 */
public class WiggleParserTest extends AbstractHeadlessTest{

    @Test
    public void testParseBedgraph_onefeat() throws Exception {
        String filepath = TestUtils.DATA_DIR + "wig/test_coords.bedgraph";
        WiggleParser parser = new WiggleParser(new ResourceLocator(filepath));
        String chr = "chr1";
        WiggleDataset ds = parser.parse();
        assertArrayEquals(new String[]{chr}, ds.getChromosomes());
        assertEquals(1.0, ds.getDataMax(), 0.01);
        assertArrayEquals(new int[]{1}, ds.getStartLocations(chr));
    }

    @Test
    public void testParseBedgraph_stats() throws Exception {

        String filepath = TestUtils.DATA_DIR + "wig/test_stats.bedgraph";
        WiggleParser parser = new WiggleParser(new ResourceLocator(filepath));

        String chr = "chr1";
        WiggleDataset ds = parser.parse();

        int len = 998;
        assertEquals(len, ds.getStartLocations(chr).length);
        int[] starts = new int[len];
        int[] ends = new int[len];
        Arrays.fill(starts, 0);
        Arrays.fill(ends, 1000);
        float[] vals = new float[]{0.15244511f, 0.9837612f, 0.17782128f, 0.81641054f, 0.89790136f};
        int[] fiStarts = ds.getStartLocations(chr);
        int[] fiEnds = ds.getEndLocations(chr);
        float[] fiVals = ds.getData("", chr);
        for(int ii=0; ii < len; ii++){
            assertEquals(starts[ii], fiStarts[ii]);
            assertEquals(ends[ii], fiEnds[ii]);
            if(ii < vals.length){
                assertEquals(vals[ii], fiVals[ii], 0.0000001f);
            }
        }

        assertEquals(0.9998707, ds.dataMax, 0.000001);

        //expected values calculated with excel which uses a different algorithm
        //Since we actually interpolate/estimate, don't expect too much precision
        assertEquals(0.086, ds.percent10, 0.001);
        assertEquals(0.89, ds.percent90, 0.01);
    }


}
