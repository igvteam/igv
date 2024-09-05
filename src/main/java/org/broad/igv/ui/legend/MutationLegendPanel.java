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

import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.PaletteColorTable;

import java.awt.*;
import java.util.Map;

/**
 * Panel to paint a legend for the mutation tracks.
 *
 * @author jrobinso
 */
public class MutationLegendPanel extends DiscreteLegendPanel {

    public MutationLegendPanel() {
        super(new PaletteColorTable());
        reloadPreferences();
    }

    @Override
    protected void reloadPreferences() {
        PaletteColorTable prefTable = PreferencesManager.getPreferences().getMutationColorScheme();
        colorTable = new PaletteColorTable();
        for (String key : prefTable.getKeys()) {
            colorTable.put(key, prefTable.get(key));
        }
        repaint();
    }

    @Override
    protected void resetPreferencesToDefault() {
        PreferencesManager.getPreferences().resetMutationColorScheme();
        reloadPreferences();
        showResetDisplay();
    }

    @Override
    protected void persistCurrentPreferences() {
        //mutation edits are persisted on edit
    }

    @Override
    public void edit() {

        boolean useColors = IGV.getInstance().getSession().getColorOverlay();
        PaletteColorTable ct = PreferencesManager.getPreferences().getMutationColorScheme();
        MutationColorMapEditor editor = new MutationColorMapEditor(IGV.getInstance().getMainFrame(), ct.getColorMap(), useColors);
        editor.setVisible(true);

        Map<String, Color> changedColors = editor.getChangedColors();
        if (!changedColors.isEmpty()) {
            for (Map.Entry<String, Color> entry : changedColors.entrySet()) {
                ct.put(entry.getKey(), entry.getValue());
            }
            String colorTableString = ct.getMapAsString();
            PreferencesManager.getPreferences().put(Constants.MUTATION_COLOR_TABLE, colorTableString);
            reloadPreferences();
        }

        boolean useColorsNew = editor.getUseColors();
        if (useColorsNew != useColors) {
            PreferencesManager.getPreferences().put(Constants.COLOR_MUTATIONS, String.valueOf(useColorsNew));
        }
    }

}
