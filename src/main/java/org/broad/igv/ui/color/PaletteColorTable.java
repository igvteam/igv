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

import org.broad.igv.Globals;

import java.awt.*;
import java.util.*;
import java.util.stream.Collectors;

/**
 * A lookup table mapping symbols (strings) -> color.  Can be initiated with our without a palette.  If
 * a palette is used colors are selected sequentially from the palette as needed until it is exhausted.
 * <p/>
 *
 * If a default color is supplied then it is used for all unknown values
 * @author jrobinso
 * @date Feb 26, 2010
 */
public class PaletteColorTable extends ColorTable {


    private static final String PALETTE_KEY = "PALETTE=";
    private static final String DEFAULT_KEY = "DEFAULT=";
    private static final Color NULL_COLOR = Color.LIGHT_GRAY;
    public static final String SERIALIZATION_ID = "PaletteColorTable";

    //list of predefined colors to use
    private final ColorPalette palette;
    private final Color defaultColor;

    public PaletteColorTable() {
        this(null, null);
    }

    public PaletteColorTable(Color defaultColor) {
        this(null, defaultColor);
    }

    public PaletteColorTable(ColorPalette palette) {
        this(palette, null);
    }

    private PaletteColorTable(ColorPalette palette, Color defaultColor){
        this.palette = palette;
        this.defaultColor = defaultColor;
    }


    public void put(String key, Color c) {
        colorMap.put(key.toLowerCase(), c);
    }

    @Override
    public Color get(String key) {
        if(key == null) {
             return NULL_COLOR;
        }
        return super.get(key.toLowerCase());
    }

    @Override
    protected Color computeColor(String key) {
        final Color c;
        if(defaultColor != null) {
            c = defaultColor;
        } else {
            final int colorIdx = colorMap.size();
            if (palette != null && colorIdx < palette.colors().length) {
                c = palette.colors()[colorIdx];
            } else {
                c = ColorUtilities.randomColor(colorIdx);
            }
        }
        return c;
    }

    @Override
    public String asString(){
        StringBuilder sb = new StringBuilder(SERIALIZATION_ID+";");
        sb.append(defaultColor == null ? "" : DEFAULT_KEY +ColorUtilities.colorToString(defaultColor) + ";");
        sb.append(palette == null ? "" : PALETTE_KEY + palette.name() + ";");
        sb.append(getMapAsString());
        return sb.toString();
    }

    public PaletteColorTable(String string){
        try {
            String[] tokens = Globals.semicolonPattern.split(string);
            int leadingTokens = 1;
            if (tokens.length >= leadingTokens + 1 && tokens[leadingTokens].startsWith(DEFAULT_KEY)) {
                String[] split = Globals.equalPattern.split(tokens[1]);
                defaultColor = ColorUtilities.stringToColor(split[1]);
                leadingTokens++;
            } else {
                defaultColor = null;
            }

            if (tokens.length >= leadingTokens + 1 && tokens[leadingTokens].startsWith(PALETTE_KEY)) {
                String[] split = Globals.equalPattern.split(tokens[1]);
                palette = ColorUtilities.getPalette(split[1]);
                leadingTokens++;
            } else {
                palette = null;
            }

            for (int i = leadingTokens; i < tokens.length; i++) {
                String[] kv = Globals.equalPattern.split(tokens[i]);
                colorMap.put(kv[0], ColorUtilities.stringToColor(kv[1]));
            }
        } catch( Exception e){
            throw new IllegalArgumentException("Could not parse palette color string: " + string, e);
        }
    }

    public Collection<String> getKeys() {
        return colorMap.keySet();
    }

    public Set<Map.Entry<String, Color>> entrySet() {
        return colorMap.entrySet();
    }

    public ColorPalette getPalette() {
        return palette;
    }

    public String getMapAsString() {
        return colorMap.entrySet().stream()
                .map(entry -> entry.getKey() + "=" + ColorUtilities.colorToString(entry.getValue()))
                .collect(Collectors.joining(";"));
    }

    public void restoreMapFromString(String string) {
        if (string == null || string.isEmpty()) return;
        colorMap.clear();
        String[] tokens = Globals.semicolonPattern.split(string);
        for (String t : tokens) {
            String[] kv = Globals.equalPattern.split(t);
            colorMap.put(kv[0], ColorUtilities.stringToColor(kv[1]));
        }
    }


    public Map<String,Color> getColorMap() {
        return colorMap;
    }
}
