package org.broad.igv.dev.affective;

import org.broad.igv.renderer.Renderer;
import org.broad.igv.track.AbstractTrack;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.color.ColorPalette;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.color.PaletteColorTable;
import org.broad.igv.ui.panel.MouseableRegion;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.awt.*;
import java.util.*;
import java.util.List;

/**
 * @author Jim Robinson
 * @date 1/22/12
 */
public class AffectiveAnnotationTrack extends AbstractTrack {

    int samplingRate = 8;      // <= TODO generalize this

    List<MouseableRegion> mouseRegions = new ArrayList<MouseableRegion>();

    Map<String, java.util.List<Annotation>> annotationMap;
    private int rowHeight = 12;
    private int topMargin = 10;

    ColorPalette pallete;
    PaletteColorTable colorTable;

    public AffectiveAnnotationTrack(String id, String name, Map<String, java.util.List<Annotation>> annotationMap) {
        super(id, name);
        this.annotationMap = annotationMap;
        colorTable  = new PaletteColorTable(ColorUtilities.getNextPalette());
    }

    public Renderer getRenderer() {
        return null;
    }

    @Override
    public int getHeight() {
        return topMargin + annotationMap.size() * rowHeight;
    }

    @Override
    public void renderName(Graphics2D g2D, Rectangle trackRectangle, Rectangle clipRect) {

        Font font = FontManager.getFont(rowHeight - 2);
        g2D.setFont(font);
        int rowTop = trackRectangle.y + topMargin;
        for (Map.Entry<String, List<Annotation>> entry : annotationMap.entrySet()) {
            g2D.drawString(entry.getKey(), 10, rowTop + rowHeight - 2);
            rowTop += rowHeight;
        }
    }

    /**
     * Render the track in the supplied rectangle.  It is the responsibility of the track to draw within the
     * bounds of the rectangle.
     *
     * @param context the render context
     * @param rect    the track bounds, relative to the enclosing DataPanel bounds.
     */
    public void render(RenderContext context, Rectangle rect) {

        mouseRegions.clear();

        String date = context.getChr();
        if (annotationMap == null || annotationMap.isEmpty()) return;

        int rowTop = rect.y + topMargin;
        int h = rowHeight - 2;

        double startTime = context.getOrigin();
        double endTime = context.getEndLocation();

        for (Map.Entry<String, List<Annotation>> entry : annotationMap.entrySet()) {

            int py = rowTop + 1;

            final double scale = context.getScale();
            final double origin = context.getOrigin();
            List<Annotation> annotations = entry.getValue();
            for (Annotation annot : annotations) {

                if (!date.equals(annot.date)) {
                    continue;
                }

                int ts = annot.startTime * samplingRate;
                int te = annot.endTime * samplingRate;

//                if (te < startTime) {
//                    continue;
//                } else if (ts > endTime) {
//                    break;
//                }

                int px = (int) ((ts - origin) / scale);
                int w = (int) Math.max(3, (te - ts) / scale);

                Color c = colorTable.get(annot.description);

                Graphics2D g = context.getGraphic2DForColor(c);
                g.fillRect(px, py, w, h);

                MouseableRegion mr = new MouseableRegion(new Rectangle(ts, py, (te-ts), h), "", annot.descriptiveText);
                mouseRegions.add(mr);
            }
            rowTop += rowHeight;
        }

    }

    /**
     * Return a value string for the tooltip window at the given location, or null to signal there is no value
     * at that location
     *
     * @param chr
     * @param position
     * @param y
     * @param frame
     * @return
     */
    @Override
    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {

        for(MouseableRegion mr : mouseRegions) {
            if(mr.containsPoint(position, y)) {
                return mr.getText();
            }
        }
        return null;
    }
}
