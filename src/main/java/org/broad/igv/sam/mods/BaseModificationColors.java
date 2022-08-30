package org.broad.igv.sam.mods;

import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sam.AlignmentTrack;

import java.awt.*;
import java.util.HashMap;
import java.util.Map;


/**
 * C Modifications
 * C m 5mC 5-Methylcytosine 27551
 * C h 5hmC 5-Hydroxymethylcytosine 76792
 * C f 5fC 5-Formylcytosine 76794
 * C c 5caC 5-Carboxylcytosine 76793
 * C C Ambiguity code; any C mod
 */

public class BaseModificationColors {

    private static Logger log = LogManager.getLogger(BaseModificationColors.class);


    static Color mColor = Color.red;
    static Color hColor = new Color(11, 132, 165);
    static Color oColor = new Color(111, 78, 129);
    static Color fColor = new Color(246, 200, 95);
    static Color cColor = new Color(157, 216, 102);
    static Color gColor = new Color(255, 160, 86);
    static Color eColor = new Color(141, 221, 208);
    static Color bColor = new Color(202, 71, 47);
    static Color genericColor = new Color(132, 178, 158);
    public static Color noModColor5MC = Color.blue;

    static HashMap<String, Color> colors = new HashMap<>();

    // 5MC overrides -- avoid colors close to blue for "C" modificatoins
    static HashMap<String, Color> colors5MC = new HashMap<>();

    static {
        colors.put("m", mColor);
        colors.put("h", hColor);
        colors.put("o", oColor);
        colors.put("f", fColor);
        colors.put("c", cColor);
        colors.put("g", gColor);
        colors.put("e", eColor);
        colors.put("b", bColor);
        colors5MC.put("h", new Color(255, 0, 255));
    }

    /**
     * Cache for modified colors
     */
    static Map<String, Color> modColorMap = new HashMap<>();
    static Map<String, Color> modColorMap5MC = new HashMap<>();


    public static Color getModColor(String modification, byte likelihood, AlignmentTrack.ColorOption colorOption) {

        // Note the pallete will always return a color, either an initially seeded one if supplied or a random color.
        Color baseColor = getBaseColor(modification, colorOption);

        int l = Byte.toUnsignedInt(likelihood);
        if (l > 255) {
            return baseColor;
        }

        String key = modification + "--" + l;
        if (colorOption == AlignmentTrack.ColorOption.BASE_MODIFICATION_5MC ||
                colorOption == AlignmentTrack.ColorOption.BASE_MODIFICATION_C) {

            if (!modColorMap5MC.containsKey(key)) {
                int alpha = Math.min(255, (int) (l * l / 64f - 4 * l + 256));    // quadratic
                if (l >= 128) {
                    modColorMap5MC.put(key, new Color(baseColor.getRed(), baseColor.getGreen(), baseColor.getBlue(), alpha));
                } else {
                    modColorMap5MC.put(key, new Color(noModColor5MC.getRed(), noModColor5MC.getGreen(), noModColor5MC.getBlue(), alpha));
                }
            }

            return modColorMap5MC.get(key);

        } else {
            if (l > 250) {
                return baseColor;
            }
            double threshold = 256 * PreferencesManager.getPreferences().getAsFloat("SAM.BASEMOD_THRESHOLD");
            if (l < threshold) {
                l = 0;
            }
            if (!modColorMap.containsKey(key)) {
                modColorMap.put(key, new Color(baseColor.getRed(), baseColor.getGreen(), baseColor.getBlue(), l));
            }
            return modColorMap.get(key);
        }
    }

    public static Color getNoModColor(byte likelihood) {

        // Note the pallete will always return a color, either an initially seeded one if supplied or a random color.
        Color baseColor = noModColor5MC;

        int l = Byte.toUnsignedInt(likelihood);
        if (l > 255) {
            return baseColor;
        }

        String key = "NOMOD--" + l;

        if (!modColorMap5MC.containsKey(key)) {
            int alpha = Math.min(255, (int) (l * l / 64f - 4 * l + 256));    // quadratic
            modColorMap5MC.put(key, new Color(baseColor.getRed(), baseColor.getGreen(), baseColor.getBlue(), alpha));
        }

        return modColorMap5MC.get(key);

    }

    private static Color getBaseColor(String modification, AlignmentTrack.ColorOption colorOption) {
        if ((colorOption == AlignmentTrack.ColorOption.BASE_MODIFICATION_5MC ||
                colorOption == AlignmentTrack.ColorOption.BASE_MODIFICATION_C) && colors5MC.containsKey(modification)) {
            return colors5MC.get(modification);
        } else if (colors.containsKey(modification)) {
            return colors.get(modification);
        } else {
            return genericColor;
        }
    }


}