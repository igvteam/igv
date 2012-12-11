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
package org.broad.igv.track;

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
        if(graphics != null) {
            graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        }

        this.panel = panel;
        this.graphicCacheByColor = new HashMap();
        this.referenceFrame = referenceFrame;
        this.visibleRect = visibleRect;
    }

    public Graphics2D getGraphic2DForColor(Color color) {

        Graphics2D g = graphicCacheByColor.get(color);
        if (g == null) {
            g = (Graphics2D) graphics.create();
            graphicCacheByColor.put(color, g);
            g.setColor(color);
            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

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
