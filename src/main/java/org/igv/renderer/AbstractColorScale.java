
package org.igv.renderer;

import org.igv.Globals;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.ui.UIConstants;
import org.igv.ui.color.ColorUtilities;

import java.awt.*;

/**
 * @author jrobinso
 */
public abstract class AbstractColorScale implements ColorScale {

    final protected static Color defaultColor = Color.BLACK;
    protected Color noDataColor = noDataColor();

    /**
     * Neutral (midpoint) color for a diverging heatmap scale.  Like "no data" this has to read as background --
     * a white midpoint on a dark panel reads as a bright band of signal, which is exactly backwards.
     */
    public static Color neutralColor() {
        return Globals.isDarkMode() ? UIConstants.getTrackPanelBackground() : Color.WHITE;
    }

    /**
     * "No data" has to read as background.  The stock default is a near-white gray, which on a dark panel paints
     * the whole extent of a heatmap track as a bright block.
     */
    public static Color noDataColor() {
        return PreferencesManager.getPreferences().getAsColor(Constants.NO_DATA_COLOR,
                ColorUtilities.shiftBrightness(UIConstants.getTrackPanelBackground(), 12));
    }

    public Color getColor(String symbol) {
        return defaultColor;
    }

    public Color getColor(float value) {
        return defaultColor;
    }

    /**
     * Method description
     *
     * @param color
     */
    public void setNoDataColor(Color color) {
        this.noDataColor = color;

    }

    public Color getNoDataColor() {
        return noDataColor;
    }
}
