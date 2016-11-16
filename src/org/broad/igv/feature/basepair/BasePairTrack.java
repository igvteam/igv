package org.broad.igv.feature.basepair;

import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.*;
import org.broad.igv.renderer.*;
import org.broad.igv.track.AbstractTrack;
import org.broad.igv.track.RenderContext;
import org.broad.igv.util.ResourceLocator;

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

    public BasePairTrack(ResourceLocator locator, String id, String name, Genome genome) {
        super(locator, id, name);
        basePairData = BasePairFileParser.loadData(locator, genome);
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