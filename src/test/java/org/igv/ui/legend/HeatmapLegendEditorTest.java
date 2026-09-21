package org.igv.ui.legend;

import org.igv.renderer.ContinuousColorScale;
import org.junit.Test;

import java.awt.Color;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

/**
 * The dialog itself needs a display, but the scale it builds does not.
 */
public class HeatmapLegendEditorTest {

    private static final Color MIN = new Color(255, 255, 204);
    private static final Color MID = Color.white;
    private static final Color MAX = new Color(202, 0, 32);

    /**
     * With the double gradient box unticked the result must be a single gradient running from the chosen
     * minimum color.  Building the three color form regardless put the midpoint color at the bottom instead.
     */
    @Test
    public void testSingleGradient() {
        ContinuousColorScale scale = HeatmapLegendEditor.buildScale(false, 0, 0, 0, 2004, MIN, MID, MAX);

        assertFalse(scale.isUseDoubleGradient());
        assertEquals(0.0, scale.getMinimum(), 1e-9);
        assertEquals(2004.0, scale.getMaximum(), 1e-9);
        assertEquals(MIN, scale.getMinColor());
        assertEquals(MAX, scale.getMaxColor());
        assertEquals("The bottom of the range is the chosen minimum color", MIN, scale.getColor(0f));
    }

    /**
     * Ticking it keeps the three color form, with the range ends ordered.
     */
    @Test
    public void testDoubleGradient() {
        ContinuousColorScale scale = HeatmapLegendEditor.buildScale(true, -2, -10, 2, 10, MIN, MID, MAX);

        assertTrue(scale.isUseDoubleGradient());
        assertEquals(-10.0, scale.getMinimum(), 1e-9);
        assertEquals(-2.0, scale.getNegStart(), 1e-9);
        assertEquals(2.0, scale.getPosStart(), 1e-9);
        assertEquals(10.0, scale.getMaximum(), 1e-9);
        assertEquals(MID, scale.getColor(0f));
    }

    /**
     * Opening a -10..10 single gradient and pressing OK without edits must save -10..10.  A single gradient's
     * posStart is fixed at max(0, minimum), so filling the visible start field from it saved 0..10.
     */
    @Test
    public void testUneditedSingleGradientRoundTrips() {
        ContinuousColorScale scale = new ContinuousColorScale(-10, 10, MIN, MAX);
        assertEquals(0.0, scale.getPosStart(), 1e-9);    // the field the old code read

        double start = HeatmapLegendEditor.visibleRangeStart(scale);
        assertEquals(-10.0, start, 1e-9);

        ContinuousColorScale saved = HeatmapLegendEditor.buildScale(false, scale.getNegStart(),
                scale.getMinimum(), start, scale.getMaximum(), MIN, MID, MAX);
        assertEquals(-10.0, saved.getMinimum(), 1e-9);
        assertEquals(10.0, saved.getMaximum(), 1e-9);
    }

    /**
     * Reversed entries are ordered rather than producing an inverted scale.
     */
    @Test
    public void testRangeEndsOrdered() {
        ContinuousColorScale scale = HeatmapLegendEditor.buildScale(false, 0, 0, 40, 0, MIN, MID, MAX);
        assertEquals(0.0, scale.getMinimum(), 1e-9);
        assertEquals(40.0, scale.getMaximum(), 1e-9);
    }
}
