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

package org.broad.igv.ui;

import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;

import javax.swing.*;
import java.awt.*;
import java.util.Hashtable;
import java.util.Set;

/**
 * @author eflakes
 */
public class FontManager {


    private static Font defaultFont;
    static Hashtable<String, Font> fontCache = new Hashtable();
    private static double scaleFactor = 1.0;


    public static Font getDefaultFont() {
        if (defaultFont == null) {
            updateDefaultFont();
        }
        return defaultFont;
    }

    static public Font getFont(int sz) {

        final PreferenceManager prefManager = PreferenceManager.getInstance();

        int size = (int) (scaleFactor * sz);

        String fontFamily = prefManager.get(PreferenceManager.DEFAULT_FONT_FAMILY);
        int attribute = prefManager.getAsInt(PreferenceManager.DEFAULT_FONT_ATTRIBUTE);
        String key = fontFamily + "_" + attribute + "_" + size;
        Font font = fontCache.get(key);
        if (font == null) {
            font = new Font(fontFamily, attribute, size);
            fontCache.put(key, font);
        }
        return font;
    }

    static public Font getFont(int attribute, int sz) {

        int size = (int) (scaleFactor * sz);

        final PreferenceManager prefManager = PreferenceManager.getInstance();
        String fontFamily = prefManager.get(PreferenceManager.DEFAULT_FONT_FAMILY);
        String key = fontFamily + "_" + attribute + "_" + size;
        Font font = fontCache.get(key);
        if (font == null) {
            font = new Font(fontFamily, attribute, size);
            fontCache.put(key, font);
        }
        return font;
    }

    public static void updateDefaultFont() {
        final PreferenceManager prefManager = PreferenceManager.getInstance();
        String fontFamily = prefManager.get(PreferenceManager.DEFAULT_FONT_FAMILY);
        int fontSize = prefManager.getAsInt(PreferenceManager.DEFAULT_FONT_SIZE);
        int attribute = prefManager.getAsInt(PreferenceManager.DEFAULT_FONT_ATTRIBUTE);
        defaultFont = new Font(fontFamily, attribute, fontSize);
    }

    public static void resetDefaultFont() {
        final PreferenceManager prefMgr = PreferenceManager.getInstance();
        prefMgr.remove(PreferenceManager.DEFAULT_FONT_SIZE);
        prefMgr.remove(PreferenceManager.DEFAULT_FONT_FAMILY);
        prefMgr.remove(PreferenceManager.DEFAULT_FONT_ATTRIBUTE);
        updateDefaultFont();
    }


    public static void updateSystemFontSize(int size) {

        Set<Object> keySet = UIManager.getLookAndFeelDefaults().keySet();
        Object[] keys = keySet.toArray(new Object[keySet.size()]);

        for (Object key : keys) {
            if (key != null && key.toString().toLowerCase().contains("font")) {
                Font font = UIManager.getDefaults().getFont(key);
                if (font != null) {
                    font = font.deriveFont((float) size);
                    UIManager.put(key, font);
                }
            }
        }
    }

    public static void scaleFontSize(double sf) {

        scaleFactor = sf;

        Set<Object> keySet = UIManager.getLookAndFeelDefaults().keySet();
        Object[] keys = keySet.toArray(new Object[keySet.size()]);
        for (Object key : keys) {
            if (key != null && key.toString().toLowerCase().contains("font")) {
                Font font = UIManager.getDefaults().getFont(key);
                if (font != null) {
                    int newSize = (int) (scaleFactor * font.getSize());
                    font = font.deriveFont((float) newSize);
                    UIManager.put(key, font);
                }
            }
        }

    }
}
