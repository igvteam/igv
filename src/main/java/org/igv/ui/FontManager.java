
package org.igv.ui;

import org.igv.prefs.IGVPreferences;
import org.igv.prefs.PreferencesManager;

import javax.swing.*;
import java.awt.*;
import java.util.Hashtable;
import java.util.Set;

import static org.igv.prefs.Constants.*;

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

        final IGVPreferences prefManager = PreferencesManager.getPreferences();

        int size = (int) (scaleFactor * sz);

        String fontFamily = prefManager.get(DEFAULT_FONT_FAMILY);
        int attribute = prefManager.getAsInt(DEFAULT_FONT_ATTRIBUTE);
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

        final IGVPreferences prefManager = PreferencesManager.getPreferences();
        String fontFamily = prefManager.get(DEFAULT_FONT_FAMILY);
        String key = fontFamily + "_" + attribute + "_" + size;
        Font font = fontCache.get(key);
        if (font == null) {
            font = new Font(fontFamily, attribute, size);
            fontCache.put(key, font);
        }
        return font;
    }

    public static void updateDefaultFont() {
        final IGVPreferences prefManager = PreferencesManager.getPreferences();
        String fontFamily = prefManager.get(DEFAULT_FONT_FAMILY);
        int fontSize = prefManager.getAsInt(DEFAULT_FONT_SIZE);
        int attribute = prefManager.getAsInt(DEFAULT_FONT_ATTRIBUTE);
        defaultFont = new Font(fontFamily, attribute, fontSize);
    }

    public static void resetDefaultFont() {
        final IGVPreferences prefMgr = PreferencesManager.getPreferences();
        prefMgr.remove(DEFAULT_FONT_SIZE);
        prefMgr.remove(DEFAULT_FONT_FAMILY);
        prefMgr.remove(DEFAULT_FONT_ATTRIBUTE);
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
