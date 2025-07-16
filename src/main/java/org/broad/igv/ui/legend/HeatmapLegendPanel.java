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

import org.broad.igv.logging.*;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import java.awt.*;
import java.text.DecimalFormat;

/**
 * @author eflakes
 */
public class HeatmapLegendPanel extends LegendPanel {

    static Logger log = LogManager.getLogger(HeatmapLegendPanel.class);


    enum Orientation {HORIZONTAL, VERTICAL}

    private Orientation orientation = Orientation.HORIZONTAL;
    private TrackType type;
    protected ContinuousColorScale colorScale;


    public HeatmapLegendPanel(TrackType type) {
        this.type = type;
        this.colorScale = PreferencesManager.getPreferences().getColorScale(type);
    }

    protected void persistResetPreferences() {
        PreferencesManager.getPreferences().setColorScale(type, colorScale);
    }

    protected void resetPreferencesToDefault() {
        colorScale = PreferencesManager.getPreferences().getDefaultColorScale(type);
        persistResetPreferences();
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

            HeatmapLegendEditor dialog = new HeatmapLegendEditor(IGV.getInstance().getMainFrame(), true, type, colorScale);

            dialog.setTitle("HeatMap Preferences");
            dialog.setVisible(true);


            if (dialog.isCanceled()) {
                IGV.getInstance().resetStatusMessage();
                return;
            }

            colorScale = dialog.getColorScheme();
            PreferencesManager.getPreferences().setColorScale(type, colorScale);
            IGV.getInstance().repaint();
            try {
                reloadPreferences();
            } finally {

                UIUtilities.invokeOnEventThread(() -> SwingUtilities.getWindowAncestor(HeatmapLegendPanel.this).toFront());
                IGV.getInstance().resetStatusMessage();
            }
        });
    }

    protected void paintLegend(Graphics2D g) {
        if (orientation == Orientation.HORIZONTAL) {
            paintHorizontal(g);
        } else {
            paintVertical(g);
        }
    }

    protected void paintHorizontal(Graphics2D g2D) {

        DecimalFormat formatter = new DecimalFormat("0.0");

        g2D.setFont(FontManager.getFont(10));

        int npts = 5;
        double max = colorScale.getMaximum();
        double min = colorScale.getMinimum();

        int w = getWidth() - 20;
        double dx = ((double) w) / npts;
        double dxj = dx / 10;
        double delta = (max - min) / npts;
        double deltaj = delta / 10;

        for (int i = 0; i < npts + 1; i++) {
            for (int j = i * 10; j < i * 10 + 10; j++) {
                double val = min + j * deltaj;

                Color c = colorScale.getColor((float) val);

                g2D.setColor(c);

                int x0 = (int) (j * dxj);
                int x1 = (int) ((j + 1) * dxj);

                g2D.fillRect(x0, 0, (x1 - x0), (int) (getHeight() / 2));
            }

            double labelVal = min + i * delta;
            int x0 = (int) (i * dx);

            g2D.setColor(Color.BLACK);
            g2D.drawString(formatter.format(labelVal), x0, (int) getHeight() - 5);
        }

    }

    void paintVertical(Graphics2D g2D) {

        DecimalFormat formatter = new DecimalFormat("0.0");

        g2D.setFont(FontManager.getFont(10));

        int npts = 5;
        double max = colorScale.getMaximum();
        double min = colorScale.getMinimum();

        int h = getWidth() - 20;
        double dy = ((double) h) / npts;
        double dyj = dy / 10;
        double delta = (max - min) / npts;
        double deltaj = delta / 10;

        int x0 = 10;
        int dx = 10;
        int y0;
        int y1 = 0;


        for (int i = 0; i < npts + 1; i++) {
            for (int j = i * 10; j < i * 10 + 10; j++) {
                double val = min + j * deltaj;

                Color c = colorScale.getColor((float) val);

                g2D.setColor(c);

                y0 = (int) (j * dyj);
                y1 = (int) ((j + 1) * dyj);

                g2D.fillRect(x0, y0, dx, y1 - y0);
            }

            double labelVal = min + i * delta;

            g2D.setColor(Color.BLACK);
            g2D.drawString(formatter.format(labelVal), x0 + 15, y1 - 5);
        }
    }

}
