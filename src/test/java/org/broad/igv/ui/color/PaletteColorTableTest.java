package org.broad.igv.ui.color;

import org.junit.Test;

import java.util.Map;

import static org.junit.Assert.assertEquals;

/**
 * @author Jim Robinson
 * @date 4/2/12
 */
public class PaletteColorTableTest {

    /**
     * Test storing a pallete color map as a string
     * @throws Exception
     */
    @Test
    public void testStoreAsString() throws Exception {

        PaletteColorTable colorTable = new PaletteColorTable(ColorUtilities.getPalette("Set 1"));
        // Force the use or creation of 50 colors
        for(int i=0; i<50; i++) {
            colorTable.get(String.valueOf(i));
        }
        String mapString = colorTable.getMapAsString();

        PaletteColorTable colorTable2 = new PaletteColorTable();
        colorTable2.restoreMapFromString(mapString);

        assertEquals(mapString, colorTable2.getMapAsString());


    }
}
