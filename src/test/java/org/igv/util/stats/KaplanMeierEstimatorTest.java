package org.igv.util.stats;

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
