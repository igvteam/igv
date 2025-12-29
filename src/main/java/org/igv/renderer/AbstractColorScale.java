/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.igv.renderer;

import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.ui.UIConstants;

import java.awt.*;

/**
 * @author jrobinso
 */
public abstract class AbstractColorScale implements ColorScale {

    final protected static Color defaultColor = Color.BLACK;
    protected Color noDataColor = PreferencesManager.getPreferences().getAsColor(Constants.NO_DATA_COLOR);

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
