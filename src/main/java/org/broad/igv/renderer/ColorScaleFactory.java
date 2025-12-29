/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.renderer;

import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 */
public class ColorScaleFactory {

    static Map<String, ColorScale> colorScaleMap = new HashMap();


    public static synchronized ColorScale getScaleFromString(String string) {

        ColorScale cs = colorScaleMap.get(string);
        if (cs == null) {

            String[] tokens = string.split(";");
            if (tokens[0].trim().equals(ContinuousColorScale.serializedClassName)) {
                cs = new ContinuousColorScale(string);
            } else if (tokens[0].trim().equals(MappedColorScale.serializationClassId)) {
                cs = new MappedColorScale(string);
            } else {
                throw new RuntimeException("Illegal ColorScale: " + string);
            }
        }
        return cs;
    }

}
