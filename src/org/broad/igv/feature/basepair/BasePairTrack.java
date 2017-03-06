package org.broad.igv.feature.basepair;

import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.*;
import org.broad.igv.renderer.*;
import org.broad.igv.session.IGVSessionReader;
import org.broad.igv.session.SubtlyImportant;
import org.broad.igv.track.AbstractTrack;
import org.broad.igv.track.RenderContext;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.collections.CollUtils;

import javax.xml.bind.annotation.*;
import java.awt.*;
import java.util.Map;

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
    private BasePairData basePairData;

    public enum ArcDirection {
        UP, DOWN
    }

    private RenderOptions renderOptions = new RenderOptions();
    
    public BasePairTrack(ResourceLocator locator, String id, String name, Genome genome) {
        super(locator, id, name);
        basePairData = BasePairFileParser.loadData(locator, genome);
    }

    public void render(RenderContext context, Rectangle rect) {
        basePairRenderer.draw(basePairData, context, rect, renderOptions);
    }

    public Renderer getRenderer() {
        return null;
    }

    public ArcDirection getArcDirection(){ return renderOptions.getArcDirection(); }

    public void setArcDirection(ArcDirection d){ renderOptions.setArcDirection(d); }

    public boolean getFitHeight() { return renderOptions.getFitHeight(); }

    public void setFitHeight(boolean b) { renderOptions.setFitHeight(b); }

    @SubtlyImportant
    private RenderOptions getRenderOptions() {
        return this.renderOptions;
    }

    @XmlElement(name = RenderOptions.NAME)
    private void setRenderOptions(RenderOptions renderOptions) {
        this.renderOptions = renderOptions;
    }

    @XmlType(name = RenderOptions.NAME)
    @XmlAccessorType(XmlAccessType.NONE)
    public static class RenderOptions {

        public static final String NAME = "BPRenderOptions";

        @XmlAttribute
        ArcDirection arcDirection;
        @XmlAttribute
        boolean fitHeight;


        RenderOptions() {
            // TODO: load some options from global PreferenceManager like AlignmentTrack?
            arcDirection = ArcDirection.DOWN;
            fitHeight = false; // scale arc heights to fit current track height
            System.out.println(">>>>>>>  BasePairTrack.RenderOptions bare constructor called <<<<<<<");
        }

        // copied these from AlignmentTrack.RenderOptions, might not be used?
        private <T extends Enum<T>> T getFromMap(Map<String, String> attributes, String key, Class<T> clazz, T defaultValue) {
            String value = attributes.get(key);
            if (value == null) {
                return defaultValue;
            }
            return CollUtils.<T>valueOf(clazz, value, defaultValue);
        }

        private String getFromMap(Map<String, String> attributes, String key, String defaultValue) {
            String value = attributes.get(key);
            if (value == null) {
                return defaultValue;
            }
            return value;
        }

        public boolean getFitHeight() { return fitHeight; }

        public void setFitHeight(boolean b) { this.fitHeight = b; }

        public ArcDirection getArcDirection() { return arcDirection; }

        public void setArcDirection(ArcDirection d) { this.arcDirection = d; }
    }

    @SubtlyImportant
    private static BasePairTrack getNextTrack() {
        return (BasePairTrack) IGVSessionReader.getNextTrack();
    }
}