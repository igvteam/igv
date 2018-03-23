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
