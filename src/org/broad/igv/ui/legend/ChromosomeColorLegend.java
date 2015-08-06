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
        paintLegend(g);
    }


    public void paintLegend(Graphics g) {

        Graphics textGraphics = g.create();
        if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.ENABLE_ANTIALISING)) {
            ((Graphics2D) textGraphics).setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        }
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
