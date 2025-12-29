/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.igv.renderer;

import org.junit.AfterClass;
import static org.junit.Assert.assertEquals;
import org.junit.BeforeClass;
import org.junit.Test;

import java.awt.*;
import java.util.HashMap;

/**
 * @author jrobinso
 */
public class MappedColorScaleTest {

    public MappedColorScaleTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    /**
     * Test of asString method, of class MappedColorScale.
     */
    @Test
    public void testSerialization() {

        HashMap<String, Color> colorMap = new HashMap();
        colorMap.put("abc", Color.black);
        colorMap.put("def", Color.blue);
        colorMap.put("ghi", Color.red);

        MappedColorScale cs = new MappedColorScale(colorMap);

        String stringRep = cs.asString();

        MappedColorScale cs2 = new MappedColorScale(stringRep);

        assertEquals(colorMap.size(), cs.getSize());
        assertEquals(cs.getSize(), cs2.getSize());

        for (String key : colorMap.keySet()) {
            assertEquals(colorMap.get(key), cs.getColor(key));
            assertEquals(cs.getColor(key), cs2.getColor(key));
        }

    }


}