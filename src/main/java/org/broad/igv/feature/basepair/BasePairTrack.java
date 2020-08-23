package org.broad.igv.feature.basepair;

import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.*;
import org.broad.igv.renderer.*;
import org.broad.igv.session.Persistable;

import org.broad.igv.track.AbstractTrack;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Show base-pairing arcs
 *
 * @author sbusan
 */

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

    public BasePairTrack() {
        this.genome = GenomeManager.getInstance().getCurrentGenome();
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


    private void setRenderOptions(BasePairTrack.RenderOptions renderOptions) {
        this.renderOptions = renderOptions;
    }

    public BasePairTrack.RenderOptions getRenderOptions() {
        return this.renderOptions;
    }

    public static class RenderOptions implements Persistable  {

        public static final String NAME = "BPRenderOptions";

        private ArcDirection arcDirection;

        private List<String> colors;

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

        @Override
        public void marshalXML(Document document, Element element) {

            if(arcDirection != ArcDirection.DOWN) {
                element.setAttribute("arcDirection", arcDirection.toString());
            }
            if(colors.size() > 0) {
                for(String c : colors) {
                    Element colorElement = document.createElement("colors");
                    colorElement.setTextContent(c);
                    element.appendChild(colorElement);
                }
            }
        }

        @Override
        public void unmarshalXML(Element element, Integer version) {

            if(element.hasAttribute("arcDirection")) {
                this.arcDirection = ArcDirection.valueOf(element.getAttribute("arcDirection"));
            }

            NodeList colorList = element.getElementsByTagName("colors");
            if(colorList.getLength() > 0) {
                this.colors = new ArrayList<>();
                for(int i=0; i<colorList.getLength(); i++) {

                    Node node = colorList.item(i);
                            this.colors.add(node.getTextContent());


                }
            }
        }
    }

    public Renderer getRenderer() {
        return null;
    }

    public void marshalXML(Document document, Element element) {

        super.marshalXML(document, element);

        Element renderElement = document.createElement("BPRenderOptions");
        renderOptions.marshalXML(document, renderElement);

        element.appendChild(renderElement);

    }

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);

        NodeList tmp = element.getElementsByTagName("BPRenderOptions");
        if(tmp.getLength() > 0) {
            Element renderElement = (Element) tmp.item(0);
            this.renderOptions.unmarshalXML(renderElement, version);
        }


    }
}