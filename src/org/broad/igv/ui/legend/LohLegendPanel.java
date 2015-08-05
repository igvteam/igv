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

import org.broad.igv.PreferenceManager;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.FontManager;

import java.awt.*;

/**
 * @author jrobinso
 */
public class LohLegendPanel extends HeatmapLegendPanel {


    public LohLegendPanel() {
        super(TrackType.LOH);
    }


    @Override
    public void paintLegend(Graphics g) {

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
            int y = getHeight() / 2;

            String label = "Loss";
            int labelWidth = (int) fm.getStringBounds(label, g2D).getWidth();
            g2D.setColor(colorScale.getMaxColor());
            g2D.fillRect(x, y, 10, 10);
            g2D.setColor(Color.BLACK);
            g2D.drawRect(x, y, 10, 10);
            g2D.drawString(label, x + 20, y + dh);
            x += labelWidth + 60;

            label = "Retained";
            labelWidth = (int) fm.getStringBounds(label, g2D).getWidth();
            g2D.setColor(colorScale.getMidColor());
            g2D.fillRect(x, y, 10, 10);
            g2D.setColor(Color.BLACK);
            g2D.drawRect(x, y, 10, 10);
            g2D.drawString(label, x + 20, y + dh);
            x += labelWidth + 60;

            label = "Conflict";
            labelWidth = (int) fm.getStringBounds(label, g2D).getWidth();
            g2D.setColor(colorScale.getMinColor());
            g2D.fillRect(x, y, 10, 10);
            g2D.setColor(Color.BLACK);
            g2D.drawRect(x, y, 10, 10);
            g2D.drawString(label, x + 20, y + dh);
            x += labelWidth + 60;

        } finally {
            g2D.dispose();
        }
    }
}
