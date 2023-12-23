package org.broad.igv.track;

import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.util.List;

public class SegTrack extends CompositeTrack {
    public SegTrack(ResourceLocator locator, List<Track> tracks) {
        super(locator, tracks);
    }

    @Override
    public void renderName(Graphics2D g2D, Rectangle trackRectangle, Rectangle visibleRectangle) {
        int trackY = 0;
        for (Track t : tracks) {
            Rectangle r = new Rectangle(trackRectangle);
            r.y = trackY;
            r.height = t.getContentHeight();
            t.renderName(g2D, r, visibleRectangle);
            trackY += r.height;
        }
    }
}
