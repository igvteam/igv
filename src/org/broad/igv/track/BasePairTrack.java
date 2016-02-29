package org.broad.igv.track;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.basepair.BasePairData;
import org.broad.igv.renderer.*;
import org.broad.igv.session.IGVSessionReader;
import org.broad.igv.session.SessionXmlAdapters;
import org.broad.igv.session.SubtlyImportant;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlType;
import javax.xml.bind.annotation.adapters.XmlJavaTypeAdapter;
import java.awt.*;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;

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