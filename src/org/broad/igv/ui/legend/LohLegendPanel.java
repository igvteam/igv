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
            g2D.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
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
