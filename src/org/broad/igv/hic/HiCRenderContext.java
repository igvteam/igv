package org.broad.igv.hic;

import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.panel.ReferenceFrame;

import javax.swing.*;
import java.awt.*;

/**
 * @author Jim Robinson
 * @date 4/13/12
 */
public class HiCRenderContext implements RenderContext {

    JComponent parent;
    Graphics2D graphics;
    Context context;
    Rectangle visibleRect;

    public HiCRenderContext(Context context, JComponent panel,
                            Graphics2D graphics, Rectangle visibleRect) {
        this.parent = panel;
        this.graphics = graphics;
        this.context = context;
        this.visibleRect = visibleRect;
    }

    public Graphics2D getGraphic2DForColor(Color color) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public Color getBackgroundColor() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String getChr() {
        return context.getChromosome().getName();
    }

    public double getOrigin() {
        return context.getOrigin();
    }

    public double getEndLocation() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public double getScale() {
        return context.getScale();
    }

    public Rectangle getVisibleRect() {
        return visibleRect;
    }

    public JComponent getPanel() {
        return parent;
    }

    public Graphics2D getGraphics() {
        return graphics;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public int getZoom() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String getGenomeId() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public ReferenceFrame getReferenceFrame() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public int bpToScreenPixel(double location) {
        return context.getScreenPosition(location);
    }

    public void dispose() {
        //To change body of implemented methods use File | Settings | File Templates.
    }
}
