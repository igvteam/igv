package org.broad.igv.hic.track;

import org.broad.igv.hic.Context;
import org.broad.igv.track.RenderContext;

import java.awt.*;

/**
 * @author jrobinso
 *         Date: 9/10/12
 *         Time: 3:15 PM
 */
public abstract class HiCTrack {

    private int height = 25;

    public int getHeight() {
        return height;
    }

    public abstract void render(Graphics2D g2d, Context context, Rectangle trackRectangle);
}
