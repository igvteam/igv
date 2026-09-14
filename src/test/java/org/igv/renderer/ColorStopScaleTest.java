package org.igv.renderer;

import org.junit.Test;

import java.awt.*;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class ColorStopScaleTest {

    private static final Color BLUE = new Color(0, 0, 200);
    private static final Color YELLOW = new Color(240, 200, 40);
    private static final Color RED = new Color(220, 30, 30);

    private final ColorStopScale scale = new ColorStopScale(List.of(
            new ColorStopScale.Stop(0.0001, RED),         // Out of order -- stops are sorted
            new ColorStopScale.Stop(0.05, BLUE),
            new ColorStopScale.Stop(0.01, YELLOW)));

    @Test
    public void testColorsAtStops() {
        assertEquals(RED, scale.getColor(0.0001));
        assertEquals(YELLOW, scale.getColor(0.01));
        assertEquals(BLUE, scale.getColor(0.05));
        assertEquals(List.of(0.0001, 0.01, 0.05), scale.getStops().stream().map(ColorStopScale.Stop::value).toList());
    }

    /**
     * Blending is on a log scale: the geometric mean of two stops is halfway between their colors.
     */
    @Test
    public void testLogBlending() {
        Color mid = scale.getColor(0.001);        // Halfway between 0.0001 and 0.01 on a log scale
        assertEquals(new Color(230, 115, 35), mid);
    }

    /**
     * Values outside the stops take the nearest stop's color; zero, negative (missing), and NaN values the lowest.
     */
    @Test
    public void testClamping() {
        assertEquals(BLUE, scale.getColor(0.5));
        assertEquals(RED, scale.getColor(0.000001));
        assertEquals(RED, scale.getColor(0));
        assertEquals(RED, scale.getColor(-1));
        assertEquals(RED, scale.getColor(Double.NaN));
    }

    @Test
    public void testStringForm() {
        String string = scale.asString();
        ColorScale parsed = ColorScaleFactory.getScaleFromString(string);
        assertTrue(parsed instanceof ColorStopScale);
        assertEquals(string, parsed.asString());
        assertEquals(scale.getColor(0.003), ((ColorStopScale) parsed).getColor(0.003));
    }

    @Test(expected = IllegalArgumentException.class)
    public void testStopsMustBePositive() {
        new ColorStopScale(List.of(new ColorStopScale.Stop(0, RED)));
    }
}
