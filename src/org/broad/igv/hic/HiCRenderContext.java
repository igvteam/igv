package org.broad.igv.hic;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.panel.ReferenceFrame;

import javax.swing.*;
import java.awt.*;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Jim Robinson
 * @date 4/13/12
 */
public class HiCRenderContext implements RenderContext {

    JComponent parent;
    Graphics2D graphics;
    Context context;
    Rectangle visibleRect;
    ReferenceFrame referenceFrame;

    private Map<Color, Graphics2D> graphicCacheByColor;

    public HiCRenderContext(Context context, JComponent panel,
                            Graphics2D graphics, Rectangle visibleRect,
                            Genome genome) {
        this.parent = panel;
        this.graphics = graphics;
        this.context = context;
        this.visibleRect = visibleRect;
        this.referenceFrame = new HiCReferenceFrame("HiC", genome);

        this.graphicCacheByColor = new HashMap();
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
        return context.getChromosomePosition(500);  //To change body of implemented methods use File | Settings | File Templates.
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
        return referenceFrame;
    }

    public int bpToScreenPixel(double location) {
        return context.getScreenPosition(location);
    }


    public Graphics2D getGraphic2DForColor(Color color) {

        Graphics2D g = graphicCacheByColor.get(color);
        if (g == null) {
            g = (Graphics2D) graphics.create();
            graphicCacheByColor.put(color, g);
            g.setColor(color);
        }
        return g;
    }

    public void dispose() {
        for (Graphics2D g : graphicCacheByColor.values()) {
            g.dispose();
        }
        graphicCacheByColor.clear();
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


}
