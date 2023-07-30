package org.broad.igv.sam.mods;

import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
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

    static HashMap<String, Color> colors = new HashMap<>();

    // 5MC overrides -- avoid colors close to blue for "C" modificatoins
    static HashMap<String, Color> colors5MC = new HashMap<>();

    static Color genericColor = new Color(132, 178, 158);
    public static Color noModColor = Color.blue;

    static Color hColor = new Color(11, 132, 165);
    static Color oColor = new Color(111, 78, 129);
    static Color fColor = new Color(246, 200, 95);
    static Color cColor = new Color(157, 216, 102);
    static Color gColor = new Color(255, 160, 86);
    static Color eColor = new Color(141, 221, 208);
    static Color bColor = new Color(0, 100, 47);
    static Color aColor = new Color(51, 0, 111);

    static {
        colors.put("m", Color.red);
        colors.put("h", hColor);
        colors.put("o", oColor);
        colors.put("f", fColor);
        colors.put("c", cColor);
        colors.put("g", gColor);
        colors.put("e", eColor);
        colors.put("b", bColor);
        colors.put("a", aColor);
        colors.put("NONE_A", noModColor);
        colors.put("NONE_C", noModColor);
        colors.put("NONE_T", noModColor);
        colors.put("NONE_G", noModColor);
        colors5MC.put("h", new Color(255, 0, 255));  // Modify h for 5mC to distinguish from blue
    }

    /**
     * Cache for modified colors
     */
    static Map<String, Color> modColorMap = new HashMap<>();


    public static Color getModColor(String modification, int l, AlignmentTrack.ColorOption colorOption) {

        // Note the pallete will always return a color, either an initially seeded one if supplied or a random color.
        Color baseColor = getBaseColor(modification, colorOption);

        if (l > 210) {
            return baseColor;
        }

        String key = modification + "--" + l;
        if (!modColorMap.containsKey(key)) {
            int alpha = Math.max(20, Math.min(255, (int) (l * l / 64f - 4 * l + 256)));
            modColorMap.put(key, new Color(baseColor.getRed(), baseColor.getGreen(), baseColor.getBlue(), alpha));
        }
        return modColorMap.get(key);

    }


    private static Color getBaseColor(String modification, AlignmentTrack.ColorOption colorOption) {
        if (colors.containsKey(modification)) {
            return colors.get(modification);
        } else {
            return genericColor;
        }
    }


}