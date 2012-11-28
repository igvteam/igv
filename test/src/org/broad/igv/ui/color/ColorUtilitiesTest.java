package org.broad.igv.ui.color;

import org.junit.Test;

import java.awt.*;

import static junit.framework.Assert.assertEquals;

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
}
