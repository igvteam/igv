package org.igv.renderer;

import org.junit.Test;

import java.awt.Color;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

public class ContinuousColorScaleTest {

    /**
     * A copy is the same scale.  The copy constructor used to force a double gradient, so "Set Heatmap Scale",
     * which edits a copy, turned a single gradient into a double one on OK even with no edits -- white at the
     * minimum and shifted colors throughout.
     */
    @Test
    public void testCopyKeepsGradientType() {

        ContinuousColorScale single = new ContinuousColorScale(0, 100, new Color(255, 255, 204), new Color(202, 0, 32));
        ContinuousColorScale copy = new ContinuousColorScale(single);

        assertFalse(copy.isUseDoubleGradient());
        assertEquals(single.asString(), copy.asString());
        for (float v : new float[]{0, 1, 25, 50, 100}) {
            assertEquals("value " + v, single.getColor(v), copy.getColor(v));
        }

        ContinuousColorScale doubleGradient =
                new ContinuousColorScale(-2, -10, 2, 10, Color.blue, Color.white, Color.red);
        ContinuousColorScale doubleCopy = new ContinuousColorScale(doubleGradient);
        assertTrue(doubleCopy.isUseDoubleGradient());
        assertEquals(doubleGradient.asString(), doubleCopy.asString());
    }
}
