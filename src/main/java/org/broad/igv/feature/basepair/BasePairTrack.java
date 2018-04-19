package org.broad.igv.feature.basepair;

import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.*;
import org.broad.igv.renderer.*;
import org.broad.igv.session.IGVSessionReader;
import org.broad.igv.session.SubtlyImportant;
import org.broad.igv.track.AbstractTrack;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.panel.ReferenceFrame;
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
    Genome genome;

    public enum ArcDirection {
        UP, DOWN
    }

    private RenderOptions renderOptions = new RenderOptions();

    public BasePairTrack(ResourceLocator locator, String id, String name, Genome genome) {
        super(locator, id, name);
        BasePairFileParser.loadData(locator, genome,
                basePairData, renderOptions);
        this.genome = genome;
    }

    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {
        return basePairData != null;
    }

    @Override
    public void load(ReferenceFrame frame) {
        BasePairFileParser.loadData(this.getResourceLocator(),
                genome, basePairData, renderOptions);
    }

    public void render(RenderContext context, Rectangle rect) {

        Graphics2D g2d = context.getGraphics();
        Rectangle clip = new Rectangle(g2d.getClip().getBounds());
        g2d.setClip(rect.intersection(clip.getBounds()));
        context.clearGraphicsCache();

        try {
            basePairRenderer.draw(basePairData, context, rect, renderOptions);
            context.clearGraphicsCache();
        } finally {
            g2d.setClip(clip);
        }
    }

    @XmlElement(name = RenderOptions.NAME)
    private void setRenderOptions(BasePairTrack.RenderOptions renderOptions) {
        this.renderOptions = renderOptions;
    }

    @SubtlyImportant
    public BasePairTrack.RenderOptions getRenderOptions() {
        return this.renderOptions;
    }

    @XmlType(name = RenderOptions.NAME)
    @XmlAccessorType(XmlAccessType.NONE)
    public static class RenderOptions {

        public static final String NAME = "BPRenderOptions";

        @XmlAttribute
        private ArcDirection arcDirection;
        // FIXME: instead of strings use @XmlJavaTypeAdapter(SessionXmlAdapters.Color.class)
        @XmlElement
        private List<String> colors;
        @XmlElement
        private List<String> colorLabels; // menu legend labels for each color

        public RenderOptions() {
            // TODO: load some options from global PreferenceManager like AlignmentTrack does?
            arcDirection = ArcDirection.DOWN;
            colors = new ArrayList();
            colorLabels = new ArrayList();
        }

        public void changeColor(Color currentColor, String currentLabel, Color newColor) {
            String currentColorString = ColorUtilities.colorToString(currentColor);
            String newColorString = ColorUtilities.colorToString(newColor);
            for (int i=0; i<getColors().size(); ++i) {
                String colorString = getColors().get(i);
                String label = getColorLabels().get(i);
                if (!currentColorString.equals(colorString) || !currentLabel.equals(label)) {
                    continue;
                }
                setColor(i, newColorString);
            }
        }

        public ArcDirection getArcDirection() {
            return arcDirection;
        }

        public void setArcDirection(ArcDirection d) {
            this.arcDirection = d;
        }

        public List<String> getColors() {
            return this.colors;
        }

        public void setColors(List<String> l) {
            this.colors = l;
        }

        public List<String> getColorLabels() {
            return this.colorLabels;
        }

        public void setColorLabels(List<String> l) {
            this.colorLabels = l;
        }

        public String getColor(int i) {
            return this.colors.get(i);
        }

        public void setColor(int i, String s) {
            this.colors.set(i, s);
        }

        public String getColorLabel(int i) {
            return this.colorLabels.get(i);
        }

        public void setColorLabel(int i, String s) {
            this.colorLabels.set(i, s);
        }

    }

    public Renderer getRenderer() {
        return null;
    }

    @SubtlyImportant
    private static BasePairTrack getNextTrack() {
        return (BasePairTrack) IGVSessionReader.getNextTrack();
    }
}