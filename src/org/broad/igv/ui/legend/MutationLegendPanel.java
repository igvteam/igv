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
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.util.ColorTable;
import org.broad.igv.ui.util.PropertyDialog;
import org.broad.igv.ui.util.PropertyDialog.PreferenceDescriptor;

import static org.broad.igv.ui.util.PropertyDialog.PreferenceType.COLOR;

import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import java.awt.*;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Panel to paint a legend for the mutation tracks.
 * TODO -- combine with other legend panels in a general
 *
 * @author jrobinso
 */
public class MutationLegendPanel extends LegendPanel {

    ColorTable colorTable;

    public MutationLegendPanel() {
        init();
    }

    private void init() {
        ColorTable prefTable = PreferenceManager.getInstance().getMutationColorScheme();
        colorTable = new ColorTable();
        for (String key : prefTable.getKeys()) {
            colorTable.put(key, prefTable.get(key));
        }
    }

    protected void persistResetPreferences() {
        PreferenceManager.getInstance().setMutationColorScheme(colorTable);
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
    public void doUserPreferences() {


        IGV.getInstance().setStatusBarMessage("Setting view properties...");

        // Add view preference items to the display
        LinkedHashMap<String, PreferenceDescriptor> labelTextToKey = addPreferences();

        Window window = SwingUtilities.getWindowAncestor(MutationLegendPanel.this);

        // Create dialog
        PropertyDialog dialog = new PropertyDialog(PreferenceManager.getInstance(),
                labelTextToKey, (Dialog) window, true);

        Component parent = IGV.getMainFrame();

        dialog.setLocationRelativeTo(parent);
        dialog.setTitle("Color Preferences");
        dialog.setVisible(true);


        if (dialog.isCanceled()) {
            IGV.getInstance().resetStatusMessage();
            return;
        }

        try {
            reloadPreferences();

        }
        finally {

            UIUtilities.invokeOnEventThread(new Runnable() {
                public void run() {
                    SwingUtilities.getWindowAncestor(MutationLegendPanel.this).toFront();
                }
            });
            IGV.getInstance().resetStatusMessage();
        }

    }

    protected LinkedHashMap<String, PreferenceDescriptor> addPreferences() {

        LinkedHashMap<String, PreferenceDescriptor> labelTextToKey = new LinkedHashMap<String,
                PreferenceDescriptor>();

        addIndelColorPreference(labelTextToKey);
        addMissenseColorPreference(labelTextToKey);
        addNonsenseColorPreference(labelTextToKey);
        addSpliceSiteColorPreference(labelTextToKey);
        addSynonymousColorPreference(labelTextToKey);
        addTargetRegionColorPreference(labelTextToKey);
        addUnknownColorPreference(labelTextToKey);

        return labelTextToKey;
    }

    protected void addNonsenseColorPreference(LinkedHashMap<String, PreferenceDescriptor> labelTextToKey) {

        Color color = colorTable.get("Nonsense");
        String labelText = "Nonsense Color: ";

        labelTextToKey.put(labelText,
                new PreferenceDescriptor(PreferenceManager.MUTATION_NONSENSE_COLOR_KEY,
                        COLOR, ColorUtilities.colorToString(color)));
    }

    protected void addIndelColorPreference(LinkedHashMap<String, PreferenceDescriptor> labelTextToKey) {

        Color color = colorTable.get("Indel");
        String labelText = "Indel Color: ";

        labelTextToKey.put(labelText,
                new PreferenceDescriptor(PreferenceManager.MUTATION_INDEL_COLOR_KEY,
                        COLOR, ColorUtilities.colorToString(color)));
    }

    protected void addTargetRegionColorPreference(LinkedHashMap<String, PreferenceDescriptor> labelTextToKey) {

        Color color = colorTable.get("Targeted_Region");
        String labelText = "Target Region Color: ";

        labelTextToKey.put(labelText,
                new PreferenceDescriptor(
                        PreferenceManager.MUTATION_TARGETED_REGION_COLOR_KEY, COLOR,
                        ColorUtilities.colorToString(color)));
    }

    protected void addMissenseColorPreference(LinkedHashMap<String, PreferenceDescriptor> labelTextToKey) {

        Color color = colorTable.get("Missense");
        String labelText = "Missense Color: ";

        labelTextToKey.put(labelText,
                new PreferenceDescriptor(PreferenceManager.MUTATION_MISSENSE_COLOR_KEY,
                        COLOR, ColorUtilities.colorToString(color)));
    }

    protected void addSpliceSiteColorPreference(LinkedHashMap<String, PreferenceDescriptor> labelTextToKey) {

        Color color = colorTable.get("Splice_site");
        String labelText = "Splice Site Color: ";

        labelTextToKey.put(
                labelText,
                new PreferenceDescriptor(
                        PreferenceManager.MUTATION_SPLICE_SITE_COLOR_KEY, COLOR,
                        ColorUtilities.colorToString(color)));
    }

    protected void addSynonymousColorPreference(LinkedHashMap<String, PreferenceDescriptor> labelTextToKey) {

        Color color = colorTable.get("Synonymous");
        String labelText = "Synonymous Color: ";

        labelTextToKey.put(
                labelText,
                new PreferenceDescriptor(
                        PreferenceManager.MUTATION_SYNONYMOUS_COLOR_KEY, COLOR,
                        ColorUtilities.colorToString(color)));
    }

    protected void addUnknownColorPreference(LinkedHashMap<String, PreferenceDescriptor> labelTextToKey) {

        Color color = colorTable.get("Unknown");
        String labelText = "Unknown Color: ";

        labelTextToKey.put(labelText,
                new PreferenceDescriptor(PreferenceManager.MUTATION_UNKNOWN_COLOR_KEY,
                        COLOR, ColorUtilities.colorToString(color)));
    }

    /**
     * Method description
     *
     * @param g
     */
    @Override
    public void paintLegend(Graphics g) {


        if (colorTable == null) {
            return;
        }

        Graphics2D g2D = null;

        try {
            g2D = (Graphics2D) g.create();

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

        }
        finally {
            g2D.dispose();
        }
    }
}
