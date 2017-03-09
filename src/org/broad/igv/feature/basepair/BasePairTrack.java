package org.broad.igv.feature.basepair;

import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.*;
import org.broad.igv.renderer.*;
import org.broad.igv.session.IGVSessionReader;
import org.broad.igv.session.SubtlyImportant;
import org.broad.igv.track.AbstractTrack;
import org.broad.igv.track.RenderContext;
import org.broad.igv.util.ResourceLocator;

import javax.xml.bind.annotation.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Show base-pairing arcs
 *
 * @author sbusan
 */
@XmlType(factoryMethod = "getNextTrack")
@XmlSeeAlso(BasePairTrack.RenderOptions.class)
public class BasePairTrack extends AbstractTrack {

    private static Logger log = Logger.getLogger(BasePairTrack.class);

    private BasePairRenderer basePairRenderer = new BasePairRenderer();
    private BasePairData basePairData = new BasePairData();

    public enum ArcDirection {
        UP, DOWN
    }

    private RenderOptions renderOptions = new RenderOptions();

    public BasePairTrack(ResourceLocator locator, String id, String name, Genome genome) {
        super(locator, id, name);
        BasePairFileParser.loadData(locator, genome,
                                    basePairData, renderOptions);
        // WIP: store colors and color labels in RenderOptions
    }

    public void render(RenderContext context, Rectangle rect) {
        basePairRenderer.draw(basePairData, context, rect, renderOptions);
    }

    public Renderer getRenderer() {
        return null;
    }

    /*public ArcDirection getArcDirection(){ return renderOptions.getArcDirection(); }

    public void setArcDirection(ArcDirection d){ renderOptions.setArcDirection(d); }

    public boolean getFitHeight() { return renderOptions.getFitHeight(); }

    public void setFitHeight(boolean b) { renderOptions.setFitHeight(b); }*/

    public RenderOptions getRenderOptions() {
        return this.renderOptions;
    }

    @XmlElement(name = RenderOptions.NAME)
    public void setRenderOptions(RenderOptions renderOptions) {
        this.renderOptions = renderOptions;
    }

    @XmlType(name = RenderOptions.NAME)
    @XmlAccessorType(XmlAccessType.NONE)
    public static class RenderOptions {

        public static final String NAME = "BPRenderOptions";

        @XmlAttribute
        ArcDirection arcDirection;
        @XmlAttribute
        boolean fitHeight; // scale arc heights to fit current track height
        @XmlAttribute
        List<String> colors; // needs to be String, not Color so XML conversion works happily
        @XmlAttribute
        List<String> colorLabels; // menu legend labels for each color

        RenderOptions() {
            // TODO: load some options from global PreferenceManager like AlignmentTrack does?
            arcDirection = ArcDirection.DOWN;
            fitHeight = false;
            colors = new ArrayList();
            colorLabels = new ArrayList();
        }

        public boolean getFitHeight() { return fitHeight; }

        public void setFitHeight(boolean b) { this.fitHeight = b; }

        public ArcDirection getArcDirection() { return arcDirection; }

        public void setArcDirection(ArcDirection d) { this.arcDirection = d; }

        public List<String> getColors() { return this.colors; }

        public void setColors(List<String> l) { this.colors = l; }

        public List<String> getColorLabels() { return this.colorLabels; }

        public void setColorLabels(List<String> l) { this.colorLabels = l; }

        public String getColor(int i) { return this.colors.get(i); }

        public void setColor(int i, String s) { this.colors.set(i, s); }

        public String getColorLabel(int i) { return this.colorLabels.get(i); }

        public void setColorLabel(int i, String s) { this.colorLabels.set(i, s); }


    }

    @SubtlyImportant
    private static BasePairTrack getNextTrack() {
        return (BasePairTrack) IGVSessionReader.getNextTrack();
    }
}