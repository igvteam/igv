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

package org.broad.igv.tools;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.sam.AlignmentCounts;
import org.broad.igv.sam.DenseAlignmentCounts;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.util.Arrays;
import java.util.Collection;

import static junit.framework.Assert.assertEquals;

/**
 * @author jacob
 * @date 2013-Jun-24
 */
@RunWith(Parameterized.class)
public class CavenerConsensusCalculatorTest extends AbstractHeadlessTest {

    private String inputBases;
    private char expectedConsensus;

    public CavenerConsensusCalculatorTest(String inputBases, char expectedConsensus) {
        this.inputBases = inputBases;
        this.expectedConsensus = expectedConsensus;

        //System.out.println("Input: " + this.inputBases);
    }

    @Test
    public void testCalcCorrectly() throws Exception {
        CavenerConsensusCalculator calc = new CavenerConsensusCalculator();
        AlignmentCounts counts = createCounts();
        assertEquals(this.expectedConsensus, calc.calculateConsensusBase(counts, 0));
    }

    private AlignmentCounts createCounts() {
        MockAlignmentCounts counts = new MockAlignmentCounts(0, 1);
        for (char c : this.inputBases.toCharArray()) {
            counts.incPositionCount(0, (byte) c);
        }
        counts.finish();

        return counts;
    }


    @Parameterized.Parameters
    public static Collection<Object[]> data() {
        Object[][] data = new Object[][]{
                {"aaa", 'a'}, {"ccccccg", 'c'}, {"g", 'g'}, {"n", 'n'},
                {"aac", 'a'}, {"acgt", 'n'},
                {"aacc", 'm'}, {"aacct", 'm'}, {"ggaaactgga", 'r'}, {"ggggccccatat", 'n'}


        };
        return Arrays.asList(data);
    }

    private static class MockAlignmentCounts extends DenseAlignmentCounts {

        public MockAlignmentCounts(int start, int end) {
            super(start, end, null);
        }

        public void incPositionCount(int pos, byte b) {
            super.incPositionCount(pos, b, (byte) 0, false);
        }
    }
}
