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

}
