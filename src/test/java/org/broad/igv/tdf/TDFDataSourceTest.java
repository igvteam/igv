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

package org.broad.igv.tdf;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.Globals;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.ChromosomeCoordinate;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static junit.framework.Assert.assertEquals;

/**
 * User: jacob
 * Date: 2012-Dec-28
 */
public class TDFDataSourceTest extends AbstractHeadlessTest {

    @Test
    public void testWholeGenomeScores() throws Exception {


        Map<String, Double> expectedValues = new HashMap<String, Double>();
        String[] chrNames = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
                "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
                "chrX", "chrY", "chrM"};
        double value = 1;
        for (String chr : chrNames) {
            expectedValues.put(chr, value);
            value++;
        }

        String testFile = TestUtils.DATA_DIR + "tdf/wholeGenomeTest.tdf";
        TDFReader reader = new TDFReader(new ResourceLocator(testFile));
        TDFDataSource ds = new TDFDataSource(reader, 0, "", genome);

        List<LocusScore> wgScores = ds.getSummaryScores(Globals.CHR_ALL, 0, Integer.MAX_VALUE, 0);
        for (LocusScore score : wgScores) {

            int genomeStart = score.getStart();
            ChromosomeCoordinate stCoord = genome.getChromosomeCoordinate(genomeStart);
            ChromosomeCoordinate endCoord = genome.getChromosomeCoordinate(score.getEnd());

            // Skip "edge" bins that include data from 2 chromosomes.  We don't know the expected value for these.
            if (stCoord.getChr().equals(endCoord.getChr())) {
                String chrName = stCoord.getChr();
                double ev = expectedValues.get(chrName);
                assertEquals(ev, score.getScore(), 0.001);
            }

        }

    }
}
