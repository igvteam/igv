package org.igv.circview.util;

import org.junit.Test;

import java.awt.Color;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

public class ColorUtilsTest {

    @Test
    public void parsesRgb() {
        Color c = ColorUtils.parseColor("rgb(80, 80, 255)");
        assertEquals(new Color(80, 80, 255), c);
        assertEquals(255, c.getAlpha());
    }

    @Test
    public void parsesRgbaWithFractionalAlpha() {
        Color c = ColorUtils.parseColor("rgba(0, 0, 255, 0.5)");
        assertEquals(0, c.getRed());
        assertEquals(0, c.getGreen());
        assertEquals(255, c.getBlue());
        assertEquals(Math.round(0.5 * 255), c.getAlpha());
    }

    @Test
    public void parsesHex() {
        assertEquals(new Color(0xAA, 0xBB, 0xCC), ColorUtils.parseColor("#aabbcc"));
        assertEquals(new Color(0xAA, 0xBB, 0xCC), ColorUtils.parseColor("#abc"));
    }

    @Test
    public void parsesNamedColor() {
        assertEquals(Color.BLACK, ColorUtils.parseColor("black"));
    }

    @Test
    public void alphaRoundTrip() {
        Color base = new Color(10, 20, 30);
        Color faded = ColorUtils.setAlpha(base, 0.25f);
        assertEquals(Math.round(0.25f * 255), faded.getAlpha());
        assertEquals(0.25f, ColorUtils.getAlpha(faded), 1f / 255f);
        // RGB preserved
        assertEquals(10, faded.getRed());
        assertEquals(20, faded.getGreen());
        assertEquals(30, faded.getBlue());
    }

    @Test
    public void rejectsGarbage() {
        assertThrows(IllegalArgumentException.class, () -> ColorUtils.parseColor("not-a-color"));
    }
}
