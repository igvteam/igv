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

package org.broad.igv.util;

import org.broad.igv.ui.color.ColorUtilities;

import java.awt.*;
import java.util.HashMap;
import java.util.Map;

public class ChromosomeColors {

    private static Map<String, Color> colorMap ;


    public static synchronized Color getColor(String chr) {
        if (!chr.startsWith("chr")) {
            chr = "chr" + chr;
        }
        if(colorMap == null) {
            initColorMap();
        }
        if (colorMap.containsKey(chr)) {
            return colorMap.get(chr);
        }
        else {
            Color color = ColorUtilities.randomColor(colorMap.size());
            colorMap.put(chr, color);
            return color;
        }
    }


    private static void initColorMap() {
        colorMap = new HashMap();
        colorMap.put("chrX", new Color(204, 153, 0));
        colorMap.put("chrY", new Color(153, 204, 0));
        colorMap.put("chrUn", Color.DARK_GRAY);
        colorMap.put("chr1", new Color(80, 80, 255));
        //colorMap.put("chr1", Color.red);
        colorMap.put("chrI", new Color(139, 155, 187));
        colorMap.put("chr2", new Color(206, 61, 50));
        colorMap.put("chrII", new Color(206, 61, 50));
        colorMap.put("chr2a", new Color(206, 61, 50).brighter());
        colorMap.put("chr2b", new Color(206, 61, 50).brighter().brighter());
        colorMap.put("chr3", new Color(116, 155, 88));
        colorMap.put("chrIII", new Color(116, 155, 88));
        colorMap.put("chr4", new Color(240, 230, 133));
        colorMap.put("chrIV", new Color(240, 230, 133));
        colorMap.put("chr5", new Color(70, 105, 131));
        colorMap.put("chr6", new Color(186, 99, 56));
        colorMap.put("chr7", new Color(93, 177, 221));
        colorMap.put("chr8", new Color(128, 34, 104));
        colorMap.put("chr9", new Color(107, 215, 107));
        colorMap.put("chr10", new Color(213, 149, 167));
        colorMap.put("chr11", new Color(146, 72, 34));
        colorMap.put("chr12", new Color(131, 123, 141));
        colorMap.put("chr13", new Color(199, 81, 39));
        colorMap.put("chr14", new Color(213, 143, 92));
        colorMap.put("chr15", new Color(122, 101, 165));
        colorMap.put("chr16", new Color(228, 175, 105));
        colorMap.put("chr17", new Color(59, 27, 83));
        colorMap.put("chr18", new Color(205, 222, 183));
        colorMap.put("chr19", new Color(97, 42, 121));
        colorMap.put("chr20", new Color(174, 31, 99));
        colorMap.put("chr21", new Color(231, 199, 111));
        colorMap.put("chr22", new Color(90, 101, 94));
        colorMap.put("chr23", new Color(204, 153, 0));
        colorMap.put("chr24", new Color(153, 204, 0));
        colorMap.put("chr25", new Color(51, 204, 0));
        colorMap.put("chr26", new Color(0, 204, 51));
        colorMap.put("chr27", new Color(0, 204, 153));
        colorMap.put("chr28", new Color(0, 153, 204));
        colorMap.put("chr29", new Color(10, 71, 255));
        colorMap.put("chr30", new Color(71, 117, 255));
        colorMap.put("chr31", new Color(255, 194, 10));
        colorMap.put("chr32", new Color(255, 209, 71));
        colorMap.put("chr33", new Color(153, 0, 51));
        colorMap.put("chr34", new Color(153, 26, 0));
        colorMap.put("chr35", new Color(153, 102, 0));
        colorMap.put("chr36", new Color(128, 153, 0));
        colorMap.put("chr37", new Color(51, 153, 0));
        colorMap.put("chr38", new Color(0, 153, 26));
        colorMap.put("chr39", new Color(0, 153, 102));
        colorMap.put("chr40", new Color(0, 128, 153));
        colorMap.put("chr41", new Color(0, 51, 153));
        colorMap.put("chr42", new Color(26, 0, 153));
        colorMap.put("chr43", new Color(102, 0, 153));
        colorMap.put("chr44", new Color(153, 0, 128));
        colorMap.put("chr45", new Color(214, 0, 71));
        colorMap.put("chr46", new Color(255, 20, 99));
        colorMap.put("chr47", new Color(0, 214, 143));
        colorMap.put("chr48", new Color(20, 255, 177));

    }
}

