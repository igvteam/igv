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

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.sam.PEStats;
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

    static Logger log = Logger.getLogger(HeatmapLegendPanel.class);


    enum Orientation {HORIZONTAL, VERTICAL}

    private Orientation orientation = Orientation.HORIZONTAL;
    private TrackType type;
    protected ContinuousColorScale colorScale;


    public HeatmapLegendPanel(TrackType type) {
        this.type = type;
        this.colorScale = PreferenceManager.getInstance().getColorScale(type);
    }

    public HeatmapLegendPanel(TrackType type, Orientation orientation) {
        this(type);
        this.orientation = orientation;
    }

    protected void persistResetPreferences() {
        PreferenceManager.getInstance().setColorScale(type, colorScale);
    }

    protected void resetPreferencesToDefault() {
        // TODO -- temporary hack.  We need some specific knowledge fo the implementation
        // in order to edit it,  but do it without a cast
        colorScale = (ContinuousColorScale) PreferenceManager.getInstance().getDefaultColorScale(type);
        persistResetPreferences();
        showResetDisplay();
    }

    protected void reloadPreferences() {
        PreferenceManager.getInstance().setColorScale(type, colorScale);
        //ColorScaleFactory.clearCache();
        repaint();
    }

    /**
     * Method description
     */
    public void edit() {

        UIUtilities.invokeOnEventThread(new Runnable() {

            public void run() {

                IGV.getInstance().setStatusBarMessage("Setting view properties...");

                HeatmapLegendEditor dialog = new HeatmapLegendEditor(IGV.getMainFrame(), true, type, colorScale);

                dialog.setTitle("HeatMap Preferences");
                dialog.setVisible(true);


                if (dialog.isCanceled()) {
                    IGV.getInstance().resetStatusMessage();
                    return;
                }
                // TODO -- temporary hack.  We need some specific knowledge fo the implementation
                // in order to edit it,  but do it without a cast

                colorScale = dialog.getColorScheme();
                PreferenceManager.getInstance().setColorScale(type, colorScale);
                IGV.getInstance().repaintDataPanels();
                try {

                    reloadPreferences();

                } finally {

                    UIUtilities.invokeOnEventThread(new Runnable() {

                        public void run() {
                            SwingUtilities.getWindowAncestor(HeatmapLegendPanel.this).toFront();
                        }
                    });
                    IGV.getInstance().resetStatusMessage();
                }
            }
        });
    }

    protected void paintLegend(Graphics g) {
        if (orientation == Orientation.HORIZONTAL) {
            paintHorizontal(g);
        } else {
            paintVertical(g);
        }
    }

    protected void paintHorizontal(Graphics g) {

        DecimalFormat formatter = new DecimalFormat("0.0");

        Graphics2D g2D = null;

        try {
            g2D = (Graphics2D) g.create();
            g2D.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
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


        } finally {
            g2D.dispose();
        }
    }

    void paintVertical(Graphics g) {
        DecimalFormat formatter = new DecimalFormat("0.0");

        Graphics2D g2D = null;

        try {
            g2D = (Graphics2D) g.create();
            g2D.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
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


        } finally {
            g2D.dispose();
        }
    }

}
