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

import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import java.awt.*;
import java.text.DecimalFormat;

/**
 * @author eflakes
 */
public class ContinuousLegendPanel extends LegendPanel {

    private final String name;

    protected ContinuousColorScale colorScale;

    public ContinuousLegendPanel(String name, ContinuousColorScale scale) {
        this.name = name;
        this.colorScale = scale;
    }

    public ContinuousColorScale getColorScale() {
        return colorScale;
    }

    protected void persistCurrentPreferences() {
        //PreferencesManager.getPreferences().setColorScale(type, colorScale);
    }

    protected void resetPreferencesToDefault() {
//        // TODO -- temporary hack.  We need some specific knowledge fo the implementation
//        // in order to edit it,  but do it without a cast
//        colorScale = IGVPreferences.getDefaultColorScale(type);
//        persistCurrentPreferences();
//        showResetDisplay();
   }

    protected void reloadPreferences() {
//        PreferencesManager.getPreferences().setColorScale(type, colorScale);
//        //ColorScaleFactory.clearCache();
        repaint();
    }


    public void updateColorScale(ContinuousColorScale scale){
        colorScale = scale;
        changeListeners.forEach(c -> c.accept(colorScale));
        repaint();
    }


    /**
     * Method description
     */
    public void edit() {

        UIUtilities.invokeOnEventThread(() -> {

            IGV.getInstance().setStatusBarMessage("Setting view properties...");

            NewContinuousLegendEditor dialog = new NewContinuousLegendEditor(IGV.getInstance().getMainFrame(), true, colorScale);

            dialog.setTitle(name + " Preferences");
            dialog.pack();
            dialog.setVisible(true);


            if (dialog.isCanceled()) {
                IGV.getInstance().resetStatusMessage();
                return;
            }

            colorScale = dialog.getColorScale();
            changeListeners.forEach(c -> c.accept(colorScale));
            //PreferencesManager.getPreferences().setColorScale(type, colorScale);
            IGV.getInstance().repaint();
            try {
                reloadPreferences();
            } finally {
                UIUtilities.invokeOnEventThread(() -> SwingUtilities.getWindowAncestor(ContinuousLegendPanel.this).toFront());
                IGV.getInstance().resetStatusMessage();
            }
        });
    }

    protected void paintLegend(Graphics2D g) {

        final DecimalFormat formatter = new DecimalFormat("0.0");

//        final Font font = FontManager.getFont(10);
//        g.setFont(font);

        final int textHeight = g.getFontMetrics().getHeight();
        final int compactHeight = textHeight * 2;
        final int npts = 5;
        final double max = colorScale.getMaximum();
        final double min = colorScale.getMinimum();

        final int w = getWidth() - 20;
        final double dx = ((double) w) / npts;
        final double dxj = dx / 10;
        final double delta = (max - min) / npts;
        final double deltaj = delta / 10;

        for (int i = 0; i < npts + 1; i++) {
            for (int j = i * 10; j < i * 10 + 10; j++) {
                double val = min + j * deltaj;

                Color c = colorScale.getColor((float) val);

                g.setColor(c);

                int x0 = (int) (j * dxj);
                int x1 = (int) ((j + 1) * dxj);

                if(getHeight() < compactHeight) { //space is limited so use all of it
                    g.fillRect(x0, 0, (x1 - x0), getHeight() );
                } else {
                    g.fillRect(x0, 0, (x1 - x0), getHeight() / 2);
                }
            }

            double labelVal = min + i * delta;
            int x0 = (int) (i * dx);
            if(getHeight() < compactHeight){ // if space is limited place text over the bar
                Color background = colorScale.getColor((float)labelVal);
                g.setColor(ColorUtilities.getVisibleOverlay(background));
                g.drawString(formatter.format(labelVal), x0, (getHeight() / 2) + (textHeight / 2));
            } else {
                g.setColor(Color.BLACK);
                g.drawString(formatter.format(labelVal), x0, getHeight() - 5);
            }
        }
    }


}
