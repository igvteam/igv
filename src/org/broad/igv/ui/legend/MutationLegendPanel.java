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

import org.broad.igv.PreferenceManager;
import org.broad.igv.renderer.ColorScale;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.PaletteColorTable;

import java.awt.*;
import java.util.Map;

/**
 * Panel to paint a legend for the mutation tracks.
 *
 * @author jrobinso
 */
public class MutationLegendPanel extends LegendPanel {

    PaletteColorTable colorTable;

    public MutationLegendPanel() {
        init();
    }

    private void init() {
        PaletteColorTable prefTable = PreferenceManager.getInstance().getMutationColorScheme();
        colorTable = new PaletteColorTable();
        for (String key : prefTable.getKeys()) {
            colorTable.put(key, prefTable.get(key));
        }
    }

    protected void persistResetPreferences() {
        PreferenceManager.getInstance().resetMutationColorScheme();
        reloadPreferences();
    }

    protected void reloadPreferences() {
        init();
        repaint();
    }

    protected ColorScale getColorScale() {

        // TODO Refactor the base class this empty method is not needed
        return null;
    }

    @Override
    protected void resetPreferencesToDefault() {

        persistResetPreferences();
        showResetDisplay();
    }

    /**
     * Open the user preferences dialog
     */
    public void edit() {

        PaletteColorTable ct = PreferenceManager.getInstance().getMutationColorScheme();
        ColorMapEditor editor = new ColorMapEditor(IGV.getMainFrame(), ct.getColorMap());
        editor.setVisible(true);

        Map<String, Color> changedColors = editor.getChangedColors();
        if (!changedColors.isEmpty()) {
            for (Map.Entry<String, Color> entry : changedColors.entrySet()) {
                ct.put(entry.getKey(), entry.getValue());
            }
            String colorTableString = ct.getMapAsString();
            PreferenceManager.getInstance().put(PreferenceManager.MUTATION_COLOR_TABLE, colorTableString);
            reloadPreferences();
        }

    }


    @Override
    public void paintLegend(Graphics g) {


        if (colorTable == null) {
            return;
        }

        Graphics2D g2D = null;

        try {
            g2D = (Graphics2D) g.create();
            if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.ENABLE_ANTIALISING)) {
                g2D.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
            }
            g2D.setFont(FontManager.getFont(10));

            FontMetrics fm = g2D.getFontMetrics();
            int dh = fm.getHeight() / 2 + 3;

            int x = 0;
            int lineHeight = 12;
            int y = lineHeight;
            int colCount = 0;

            for (Map.Entry<String, Color> entry : colorTable.entrySet()) {

                String mutType = entry.getKey();
                String label = mutType.replace("_", " ");
                int labelWidth = (int) fm.getStringBounds(label, g2D).getWidth();

                g2D.setColor(entry.getValue());
                g2D.fillRect(x, y, 10, 10);
                g2D.setColor(Color.BLACK);
                g2D.drawRect(x, y, 10, 10);
                g2D.drawString(label, x + 20, y + dh);
                x += labelWidth + 40;
                colCount++;

                if (colCount % 5 == 0) {
                    y += lineHeight + 5;
                    x = 0;
                }
            }

        } finally {
            g2D.dispose();
        }
    }
}
