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
package org.broad.igv.track;

import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.panel.ReferenceFrame;

import javax.swing.*;
import java.awt.*;
import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 */
public class RenderContextImpl implements RenderContext {

    private Graphics2D graphics;
    private Map<Color, Graphics2D> graphicCacheByColor;
    private ReferenceFrame referenceFrame;
    private JComponent panel;
    private Rectangle visibleRect;


    public RenderContextImpl(JComponent panel, Graphics2D graphics, ReferenceFrame referenceFrame, Rectangle visibleRect) {
        this.graphics = graphics;
        this.panel = panel;
        this.graphicCacheByColor = new HashMap();
        this.referenceFrame = referenceFrame;
        this.visibleRect = visibleRect;
    }

    public Graphics2D getGraphic2DForColor(Color color) {

        Graphics2D g = graphicCacheByColor.get(color);
        if (g == null) {
            g = (Graphics2D) graphics.create();
            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            graphicCacheByColor.put(color, g);
            g.setColor(color);

        }
        return g;
    }


    public Color getBackgroundColor() {
        return panel.getBackground();
    }

    public String getChr() {
        return referenceFrame.getChrName();
    }

    public double getOrigin() {
        return referenceFrame.getOrigin();
    }

    public double getEndLocation() {
        return referenceFrame.getEnd();
    }

    public double getScale() {
        return referenceFrame.getScale();
    }

    public Rectangle getVisibleRect() {
        return visibleRect;
    }

    public JComponent getPanel() {
        return panel;
    }

    public Graphics2D getGraphics() {
        return graphics;
    }

    public int getZoom() {
        return referenceFrame.getZoom();
    }

    public ReferenceFrame getReferenceFrame() {
        return referenceFrame;
    }

    public int bpToScreenPixel(double location) {
        final double scale = getScale();
        final double origin = getOrigin();
        return (int) ((location - origin) / scale);

    }

    /**
     * Release graphics objects
     *
     * @throws java.lang.Throwable
     */
    @Override
    protected void finalize() throws Throwable {
        super.finalize();
        dispose();
    }

    public void dispose() {
        for (Graphics2D g : graphicCacheByColor.values()) {
            g.dispose();
        }
        graphicCacheByColor.clear();
    }

}
