package org.broad.igv.charts;

import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * @author Jim Robinson
 * @date 10/29/11
 */
public class AxisTest {
    @Test
    public void testPerformAutoScale() throws Exception {

        double[] ticks = Axis.computeTicks(0, 0.09, 10);
        assertEquals("min tick", 0, ticks[0], 1.0e-6);
        assertEquals("min tick", 0.01, ticks[1], 1.0e-6);

        ticks = Axis.computeTicks(-0.048, 0.045, 10);
        assertEquals("min tick", -0.05, ticks[0], 1.0e-6);
        assertEquals("min tick", 0.01, ticks[1], 1.0e-6);

        ticks = Axis.computeTicks(10, 999, 10);
        assertEquals("min tick", 0, ticks[0], 1.0e-6);
        assertEquals("min tick", 100, ticks[1], 1.0e-6);


    }
}
