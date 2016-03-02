package org.broad.igv.feature.basepair;

import org.apache.log4j.Logger;
import org.broad.igv.renderer.*;
import org.broad.igv.track.AbstractTrack;
import org.broad.igv.track.RenderContext;

import java.awt.*;

/**
 * Show base-pairing arcs
 *
 * @author sbusan
 */
public class BasePairTrack extends AbstractTrack {

    private static Logger log = Logger.getLogger(BasePairTrack.class);

    private BasePairRenderer basePairRenderer = new BasePairRenderer();
    private BasePairData basePairData;

    public BasePairTrack(BasePairData data, String name) {
        super(name);
        basePairData = data;
    }

    public void render(RenderContext context, Rectangle rect) {
        basePairRenderer.draw(basePairData, context, rect);
    }

    public Renderer getRenderer() {
        return null;
    }

    public int getDirection(){
        return basePairRenderer.getDirection();
    }

    public void setDirection(int d){
        basePairRenderer.setDirection(d);
    }
}