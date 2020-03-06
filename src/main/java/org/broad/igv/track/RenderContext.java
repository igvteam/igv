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

import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sam.InsertionMarker;
import org.broad.igv.ui.panel.ReferenceFrame;

import javax.swing.*;
import java.awt.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.broad.igv.prefs.Constants.ENABLE_ANTIALISING;

/**
 * @author jrobinso
 */

public class RenderContext {

    private Graphics2D graphics;
    private Map<Object, Graphics2D> graphicCache;
    private ReferenceFrame referenceFrame;
    private JComponent panel;
    private Rectangle visibleRect;
    private boolean merged = false;
    public transient int translateX;
    private List<InsertionMarker> insertionMarkers;
    public boolean multiframe = false;


    public RenderContext(JComponent panel, Graphics2D graphics, ReferenceFrame referenceFrame, Rectangle visibleRect) {
        this.graphics = graphics;
        this.panel = panel;
        this.graphicCache = new HashMap();
        this.referenceFrame = referenceFrame;
        this.visibleRect = visibleRect;
        if (PreferencesManager.getPreferences().getAntiAliasing() && graphics != null) {
            graphics.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        }
    }

    public RenderContext(RenderContext context) {
        this.graphics = context.graphics;
        this.graphicCache = new HashMap<>();
        this.referenceFrame = new ReferenceFrame(context.referenceFrame);
        this.panel = context.panel;
        this.visibleRect = new Rectangle(context.visibleRect);
        if (PreferencesManager.getPreferences().getAntiAliasing() && graphics != null) {
            graphics.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        }
    }

    public Graphics2D getGraphics() {
        return graphics;
    }

    public void clearGraphicsCache() {
        for(Graphics2D g: graphicCache.values()) {
            g.dispose();
        }
        graphicCache.clear();
    }

    public Graphics2D getGraphics2D(Object key) {
        Graphics2D g = graphicCache.get(key);
        if (g == null) {
            g = (Graphics2D) graphics.create();
            if (PreferencesManager.getPreferences().getAntiAliasing()) {
                g.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
            }
            graphicCache.put(key, g);
        }
        return g;
    }

    public Graphics2D getGraphic2DForColor(Color color) {
        Graphics2D g = getGraphics2D(color);
        g.setColor(color);
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

    public int getZoom() {
        return referenceFrame.getZoom();
    }

    public ReferenceFrame getReferenceFrame() {
        return referenceFrame;
    }

    public void dispose() {
        // Note: don't dispose of "this.graphics", it is managed by the framwork.
        for (Graphics2D g : graphicCache.values()) {
            g.dispose();
        }
        graphicCache.clear();
    }


    public void setInsertionMarkers(List<InsertionMarker> insertionMarkers) {
        this.insertionMarkers = insertionMarkers;
    }

    public List<InsertionMarker> getInsertionMarkers() {
        return insertionMarkers;
    }
}
