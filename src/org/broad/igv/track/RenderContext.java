package org.broad.igv.track;

import org.broad.igv.ui.panel.ReferenceFrame;

import javax.swing.*;
import java.awt.*;

/**
 * @author Jim Robinson
 * @date 4/13/12
 */
public interface RenderContext {
    Graphics2D getGraphic2DForColor(Color color);

    Color getBackgroundColor();

    String getChr();

    double getOrigin();

    double getEndLocation();

    double getScale();

    Rectangle getVisibleRect();

    JComponent getPanel();

    Graphics2D getGraphics();

    int getZoom();

    ReferenceFrame getReferenceFrame();

    int bpToScreenPixel(double location);

    void dispose();

    /**
     * Data managers can load a certain number
     * of intervals at once. This sets that size.
     * Implementations may have their own floors,
     * so this may not have an effect
     */
    public void setCacheSize(int cacheSize);

    public int getCacheSize();
}
