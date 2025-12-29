package org.igv.ui.color;

import org.junit.Test;

import java.awt.*;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertNotNull;

/**
 * @author jrobinso
 *         Date: 11/28/12
 *         Time: 11:49 AM
 */
public class ColorUtilitiesTest {


    @Test
    public void testStringToColor() throws Exception {
        // Test parsing strings with quotes (common problem with Excel exports)
        String quotedString = "\"0,0,255,\"";
        Color b = ColorUtilities.stringToColor(quotedString);
        assertEquals(Color.blue, b);
    }

    @Test
    public void testHexColorString() throws Exception {
        Color grey = new Color(170,170, 170);
        Color c = ColorUtilities.stringToColor("#AAAAAA");
        assertEquals(grey, c);
    }

    @Test
    public void testRGBToColor() throws Exception {
        // Test parsing javascript style rgb string
        String rgbString = "rgb(0,0,255)";
        Color b = ColorUtilities.stringToColor(rgbString);
        assertEquals(Color.blue, b);
    }
}
