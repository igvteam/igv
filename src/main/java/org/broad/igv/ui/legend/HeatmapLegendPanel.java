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
package org.broad.igv.ui.legend;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;

/**
 * @author eflakes
 */
public class HeatmapLegendPanel extends ContinuousLegendPanel {
    private final TrackType type;

    public HeatmapLegendPanel(TrackType type) {
        super("Heatmap", PreferencesManager.getPreferences().getColorScale(type));
        this.type = type;
    }

    protected void persistCurrentPreferences() {
        PreferencesManager.getPreferences().setColorScale(type, colorScale);
    }

    protected void resetPreferencesToDefault() {
        // TODO -- temporary hack.  We need some specific knowledge fo the implementation
        // in order to edit it,  but do it without a cast
        colorScale = IGVPreferences.getDefaultColorScale(type);
        persistCurrentPreferences();
        showResetDisplay();
    }

    protected void reloadPreferences() {
        PreferencesManager.getPreferences().setColorScale(type, colorScale);
        //ColorScaleFactory.clearCache();
        repaint();
    }

    /**
     * Method description
     */
    public void edit() {

        UIUtilities.invokeOnEventThread(() -> {

            IGV.getInstance().setStatusBarMessage("Setting view properties...");

            ContinuousLegendEditor dialog = new ContinuousLegendEditor(IGV.getInstance().getMainFrame(), true, colorScale);

            dialog.setTitle("Heatmap Preferences");
            dialog.setVisible(true);


            if (dialog.isCanceled()) {
                IGV.getInstance().resetStatusMessage();
                return;
            }

            colorScale = dialog.getColorScheme();
            changeListeners.forEach(c -> c.accept(colorScale));
            //PreferencesManager.getPreferences().setColorScale(type, colorScale);
            IGV.getInstance().repaint();
            try {
                reloadPreferences();
            } finally {
                UIUtilities.invokeOnEventThread(() -> SwingUtilities.getWindowAncestor(this).toFront());
                IGV.getInstance().resetStatusMessage();
            }
        });
    }

}
