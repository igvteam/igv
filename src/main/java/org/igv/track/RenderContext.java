package org.igv.track;

import org.igv.prefs.PreferencesManager;
import org.igv.sam.InsertionMarker;
import org.igv.ui.panel.ReferenceFrame;

import javax.swing.*;
import java.awt.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.igv.prefs.Constants.ENABLE_ANTIALISING;

/**
 * @author jrobinso
 */

public class RenderContext {

    private Graphics2D graphics;
    private Map<Object, Graphics2D> graphicCache;
    private ReferenceFrame referenceFrame;
    private JComponent panel;
    public Rectangle trackRectangle;
    private Rectangle clipBounds;
    public boolean multiframe = false;
    public int expandedInsertionPosition = -1;

    /**
     * X translation for this context relative to its parent.  This is used in expanded insertion "multi-frame* view
     * to convert screen coordinates to parent reference system when recording the pixel location of drawn objects
     */
    public int translateX = 0;

    public RenderContext(JComponent panel, Graphics2D graphics, ReferenceFrame referenceFrame, Rectangle trackRectangle, Rectangle clipBounds) {
        this.graphics = graphics;
        this.panel = panel;
        this.graphicCache = new HashMap();
        this.referenceFrame = referenceFrame;
        this.trackRectangle = trackRectangle;
        this.clipBounds = clipBounds;
        if (PreferencesManager.getPreferences().getAntiAliasing() && graphics != null) {
            graphics.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        }
    }

    public RenderContext(RenderContext context) {
        this.graphics = (Graphics2D) context.graphics.create();
        this.graphicCache = new HashMap<>();
        this.referenceFrame = new ReferenceFrame(context.referenceFrame);
        this.panel = context.panel;
        this.trackRectangle = new Rectangle(context.trackRectangle);
        this.clipBounds = new Rectangle(context.clipBounds);
        this.expandedInsertionPosition = context.expandedInsertionPosition;
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
        return trackRectangle;
    }

    public Rectangle getTrackRectangle() {
        return trackRectangle;
    }

    public Rectangle getClipBounds() {
        return clipBounds;
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

}
