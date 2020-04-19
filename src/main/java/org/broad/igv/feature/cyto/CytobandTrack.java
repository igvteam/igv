package org.broad.igv.feature.cyto;

import org.broad.igv.track.AbstractTrack;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.awt.*;

public class CytobandTrack extends AbstractTrack {

    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {
        return false;
    }

    @Override
    public void load(ReferenceFrame frame) {

    }

    @Override
    public void render(RenderContext context, Rectangle rect) {

    }
}
