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
