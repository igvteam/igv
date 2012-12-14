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

import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.ui.FontManager;
import org.broad.igv.util.ChromosomeColors;

import javax.swing.*;
import java.awt.*;

/**
 * @author jrobinso
 */

public class ChromosomeColorLegend extends JPanel {


    public ChromosomeColorLegend() {
        this.setSize(480, 24);
    }

    @Override
    public Dimension getPreferredSize() {
        return new Dimension(480, 24);
    }


    /**
     * Method description
     *
     * @param g
     */
    @Override
    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        paintLegend(g);
    }


    public void paintLegend(Graphics g) {

        Graphics textGraphics = g.create();
        ((Graphics2D) textGraphics).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        textGraphics.setColor(Color.black);
        textGraphics.setFont(FontManager.getFont(10));

        int w = (int) (getWidth() / 24);
        int h = getHeight() / 2;
        int x = 0;

        for (int i = 1; i <= 24; i++) {
            String chr = (i < 23 ? "chr" + i : i == 23 ? "chrX" : "chrY");
            Color c = ChromosomeColors.getColor(chr);
            g.setColor(c);
            g.fillRect(x, 0, w, h);

            String tmp = (i < 23 ? String.valueOf(i) : i == 23 ? "X" : "Y");
            GraphicUtils.drawCenteredText(tmp, x, h, w, h, textGraphics);
            x += w;

        }


    }
}
