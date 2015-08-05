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

package org.broad.igv.ui.util;

import javax.swing.*;
import javax.swing.plaf.basic.BasicLabelUI;
import java.awt.*;
import java.awt.geom.AffineTransform;

/**
 * @author jrobinso
 * @date Sep 16, 2010
 */
public class SwitchingLabelUI extends BasicLabelUI {

    enum Orientation {
        Horizontal, Vertical
    }


    private int verticalInset;
    private Orientation orientation = null;
    private boolean clockwise = false;
    private static Rectangle paintIconR = new Rectangle();
    private static Rectangle paintTextR = new Rectangle();
    private static Rectangle paintViewR = new Rectangle();
    private static Insets paintViewInsets = new Insets(0, 0, 0, 0);

    public SwitchingLabelUI(int verticalInset) {
        this.verticalInset = verticalInset;
    }


    public Dimension getPreferredSize(JComponent c) {
        Dimension dim = super.getPreferredSize(c);
        return new Dimension(dim.height, dim.width);
    }


    public void paint(Graphics g, JComponent c) {

        JLabel label = (JLabel) c;
        if (label.getWidth() > label.getHeight()) {
            if (orientation != Orientation.Horizontal) {
                label.setHorizontalAlignment(JLabel.CENTER);
                label.setBorder(BorderFactory.createEmptyBorder(0, 0, 0, 0));
                orientation = Orientation.Horizontal;
            }
            super.paint(g, c);
        } else {
            if (orientation != Orientation.Vertical) {
                label.setHorizontalAlignment(JLabel.LEFT);
                label.setBorder(BorderFactory.createEmptyBorder(0, verticalInset, 0, 0));
                orientation = Orientation.Vertical;
            }

            String text = label.getText();
            Icon icon = (label.isEnabled()) ? label.getIcon() : label.getDisabledIcon();

            if ((icon == null) && (text == null)) {
                return;
            }

            FontMetrics fm = g.getFontMetrics();
            paintViewInsets = c.getInsets(paintViewInsets);

            paintViewR.x = paintViewInsets.left;
            paintViewR.y = paintViewInsets.top;

            // Use inverted height & width
            paintViewR.height = c.getWidth() - (paintViewInsets.left + paintViewInsets.right);
            paintViewR.width = c.getHeight() - (paintViewInsets.top + paintViewInsets.bottom);

            paintIconR.x = paintIconR.y = paintIconR.width = paintIconR.height = 0;
            paintTextR.x = paintTextR.y = paintTextR.width = paintTextR.height = 0;

            String clippedText =
                    layoutCL(label, fm, text, icon, paintViewR, paintIconR, paintTextR);

            Graphics2D g2 = (Graphics2D) g;
            AffineTransform tr = g2.getTransform();
            if (clockwise) {
                g2.rotate(Math.PI / 2);
                g2.translate(0, -c.getWidth());
            } else {
                g2.rotate(-Math.PI / 2);
                g2.translate(-c.getHeight(), 0);
            }

            if (icon != null) {
                icon.paintIcon(c, g, paintIconR.x, paintIconR.y);
            }

            if (text != null) {
                int textX = paintTextR.x;
                int textY = paintTextR.y + fm.getAscent();

                if (label.isEnabled()) {
                    paintEnabledText(label, g, clippedText, textX, textY);
                } else {
                    paintDisabledText(label, g, clippedText, textX, textY);
                }
            }


            g2.setTransform(tr);
        }
    }
}


