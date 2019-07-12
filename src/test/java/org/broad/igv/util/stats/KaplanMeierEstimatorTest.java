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

package org.broad.igv.util.stats;

import org.junit.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.*;

/**
 * @author jrobinso
 * @date Dec 4, 2010
 */
public class KaplanMeierEstimatorTest {


    /**
     * http://www.mas.ncl.ac.uk/~nmf16/teaching/mas3311/week06.pdf
     * [10 13* 18* 19 23* 30 36 38* 54* 56* 59 75 93 97 104* 107 107* 107*]
     */
    @Test
    public void testCensured() {
        int[] time = {10, 13, 18, 19, 23, 30, 36, 38, 54, 56, 59, 75, 93, 97, 104, 107, 107, 107};
        boolean[] censured = {false, true, true, false, true, false, false, true, true, true, false,
                false, false, false, true, false, true, true};
        double [] expectedSurvivals = {1.00, 0.94, 0.88, 0.81, 0.75, 0.65, 0.56, 0.47, 0.37};

        List<KaplanMeierEstimator.Interval> intervals = KaplanMeierEstimator.compute(time, censured);
        for (int j = 0; j < intervals.size(); j++) {
            KaplanMeierEstimator.Interval interval = intervals.get(j);
            assertEquals(expectedSurvivals[j], interval.getCumulativeSurvival(), 0.01);
        }

    }


   


}
