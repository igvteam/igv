/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
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
