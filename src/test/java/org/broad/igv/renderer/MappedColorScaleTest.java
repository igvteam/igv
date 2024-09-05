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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.renderer;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import java.awt.*;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.SequencedMap;

/**
 * @author jrobinso
 */
public class MappedColorScaleTest {

    @Test
    public void testNoChangesToSerialization(){
        SequencedMap<String, Color> map = new LinkedHashMap<>();
        map.put("abc", Color.black);
        map.put("def", Color.blue);
        map.put("ghi", Color.red);
        MappedColorScale expected = new MappedColorScale(map);
        String serialized = "MappedColorScale;abc 0,0,0;def 0,0,255;ghi 255,0,0";

        assertEquals(expected.asString(), serialized);
        assertColorMapsEqual(map, new MappedColorScale(serialized), expected);
    }

    /**
     * Test of asString method, of class MappedColorScale.
     */
    @Test
    public void testSerialization() {

        HashMap<String, Color> colorMap = new HashMap<>();
        colorMap.put("abc", Color.black);
        colorMap.put("def", Color.blue);
        colorMap.put("ghi", Color.red);

        MappedColorScale cs = new MappedColorScale(colorMap);

        String stringRep = cs.asString();

        MappedColorScale cs2 = new MappedColorScale(stringRep);

        assertColorMapsEqual(colorMap, cs, cs2);

    }

    private static void assertColorMapsEqual(Map<String, Color> colorMap, MappedColorScale cs, MappedColorScale cs2) {
        assertEquals(colorMap.size(), cs.getSize());
        assertEquals(cs.getSize(), cs2.getSize());

        for (String key : colorMap.keySet()) {
            assertEquals(colorMap.get(key), cs.getColor(key));
            assertEquals(cs.getColor(key), cs2.getColor(key));
        }
    }


}