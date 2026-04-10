package org.igv.feature.mut;

import org.igv.logging.*;
import org.igv.prefs.PreferencesManager;
import org.igv.renderer.Renderer;
import org.igv.track.RenderContext;
import org.igv.track.Track;
import org.igv.ui.IGV;
import org.igv.ui.FontManager;
import org.igv.ui.color.ColorPalette;
import org.igv.ui.color.ColorTable;
import org.igv.ui.color.PaletteColorTable;

import java.awt.*;
import java.util.List;

public class MutationRenderer implements Renderer<Mutation> {

    private static Logger log = LogManager.getLogger(MutationRenderer.class);

    public String getDisplayName() {
        return "Mutation";
    }

    ColorTable colorTable = PreferencesManager.getPreferences().getMutationColorScheme();

    /**
     * Note:  assumption is that featureList is sorted by start position.
     */
    public void render(List<Mutation> featureList, RenderContext context, Rectangle trackRectangle, Track track) {

        double origin = context.getOrigin();
        double locScale = context.getScale();

        Graphics2D graphics = (Graphics2D) context.getGraphics().create();

        if (featureList != null && featureList.size() > 0) {

            Rectangle lastRect = null;
            for (Mutation feature : featureList) {
                // Note -- don't cast these to an int until the range is checked.
                // could get an overflow.
                double pixelStart = ((feature.getStart() - origin) / locScale);
                double pixelEnd = ((feature.getEnd() - origin) / locScale);

                // If the any part of the feature fits in the
                // Track rectangle draw it
                if (pixelEnd >= trackRectangle.getX() && pixelStart <= trackRectangle.getMaxX()) {

                    // Set color used to draw the feature

                    Color color = feature.getColor();
                    Graphics2D g = context.getGraphic2DForColor(color);
                    g.setFont(FontManager.getDefaultFont());


                    int w = (int) (pixelEnd - pixelStart);
                    if (w < 3) {
                        w = 3;
                        pixelStart--;
                    }

                    int mutHeight = (int) Math.max(1, trackRectangle.getHeight() - 2);
                    int mutY = (int) (trackRectangle.getY() + (trackRectangle.getHeight() - mutHeight) / 2);


                    Rectangle mutRect = new Rectangle((int) pixelStart, mutY, w, mutHeight);

                    graphics.setColor(colorTable.get(feature.getMutationType()));
                    mutRect.x--;
                    mutRect.width += 2;
                    graphics.fill(mutRect);

                    lastRect = mutRect;
                }
            }

        }
        graphics.dispose();
    }
}
