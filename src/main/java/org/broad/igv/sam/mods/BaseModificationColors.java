package org.broad.igv.sam.mods;

import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sam.AlignmentRenderer;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.ui.color.ColorUtilities;

import static org.broad.igv.prefs.Constants.*;

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


    public static void updateColors() {
        IGVPreferences preferences = PreferencesManager.getPreferences();
        colors.put("m", preferences.getAsColor(BASEMOD_M_COLOR));
        colors.put("h", preferences.getAsColor(BASEMOD_H_COLOR));
        colors.put("o", preferences.getAsColor(BASEMOD_O_COLOR));
        colors.put("f", preferences.getAsColor(BASEMOD_F_COLOR));
        colors.put("c", preferences.getAsColor(BASEMOD_C_COLOR));
        colors.put("g", preferences.getAsColor(BASEMOD_G_COLOR));
        colors.put("e", preferences.getAsColor(BASEMOD_E_COLOR));
        colors.put("b", preferences.getAsColor(BASEMOD_B_COLOR));
        colors.put("a", preferences.getAsColor(BASEMOD_A_COLOR));
        colors.put("other", preferences.getAsColor(BASEMOD_OTHER_COLOR));
        colors.put("NONE_A", preferences.getAsColor(BASEMOD_NONE_A_COLOR));
        colors.put("NONE_C", preferences.getAsColor(BASEMOD_NONE_C_COLOR));
        colors.put("NONE_U", preferences.getAsColor(BASEMOD_NONE_C_COLOR));
        colors.put("NONE_T", preferences.getAsColor(BASEMOD_NONE_T_COLOR));
        colors.put("NONE_G", preferences.getAsColor(BASEMOD_NONE_G_COLOR));
        colors.put("NONE_N", preferences.getAsColor(BASEMOD_NONE_N_COLOR));
        modColorMap.clear();
    }

    /**
     * Cache for modified colors
     */
    static Map<String, Color> modColorMap = new HashMap<>();


    public static Color getModColor(String modification, int l, AlignmentTrack.ColorOption colorOption) {

        if (colors.isEmpty()) updateColors();

        // Note the pallete will always return a color, either an initially seeded one if supplied or a random color.
        Color baseColor = getBaseColor(modification);

        String key = modification + l + colorOption;
        if (!modColorMap.containsKey(key)) {
            int alpha = colorOption == AlignmentTrack.ColorOption.BASE_MODIFICATION_2COLOR ?
                    Math.max(20, Math.min(255, 20 + (int) (l * l / 50f - 4 * l + 200))) :
                    Math.max(20, (int) Math.min(255, 6.127e-3*l*l));

            modColorMap.put(key, new Color(baseColor.getRed(), baseColor.getGreen(), baseColor.getBlue(), alpha));
        }
        return modColorMap.get(key);

    }


    private static Color getBaseColor(String modification) {
        if (colors.containsKey(modification)) {
            return colors.get(modification);
        } else {
            return colors.get("other");
        }
    }


}