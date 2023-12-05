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

package org.broad.igv.feature.genome;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.load.ChromSizesParser;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

/**
 * @author jrobinso
 * Date: 2/13/13
 * Time: 2:21 PM
 */
public class ChromosomeComparatorTest {


    // Test chr names which resolve to numbers too large for "int".  In response to bug report.

    @Test
    public void testCompareLargeNumbers() throws Exception {

        String f = TestUtils.DATA_DIR + "genomes/hg19.chrom.sizes";
        List<Chromosome> chromosomeList = ChromSizesParser.parse(f);
        List<Chromosome> sortedList = new ArrayList<>(chromosomeList);

        ChromosomeComparator comp = new ChromosomeComparator();
        sortedList.sort(comp);

        for(int i=0; i<25; i++) {
            assertEquals(chromosomeList.get(i).getName(), sortedList.get(i).getName());
        }
    }
}
