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
            g2D.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
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
