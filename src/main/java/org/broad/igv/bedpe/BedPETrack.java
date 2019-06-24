package org.broad.igv.bedpe;

import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.AbstractTrack;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.TrackClickEvent;
import org.broad.igv.track.TrackMenuUtils;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ResourceLocator;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import javax.swing.*;
import java.awt.*;
import java.util.*;
import java.util.List;

import static org.broad.igv.bedpe.BedPETrack.Direction.UP;
import static org.broad.igv.track.TrackMenuUtils.refresh;


/**
 * Created by jrobinso on 6/29/18.
 */
public class BedPETrack extends AbstractTrack {

    enum Direction {UP, DOWN}

    enum GraphType {BLOCK, ARC, PROPORTIONAL_ARC}

    private Genome genome;
    BedPETrack.Direction direction = UP; //DOWN;
    GraphType graphType = GraphType.ARC;  // GraphType.block; //
    int thickness = 1;
    boolean autoscale = true;
    int gap = 5;

    private Map<String, List<BedPE>> featureMap;
    private Map<GraphType, BedPERenderer> renderers;

    public BedPETrack() {
    }

    public BedPETrack(ResourceLocator locator, List<BedPEFeature> featureList, Genome genome) {
        super(locator);
        init(featureList, genome);
        this.genome = genome;
        setHeight(250, true);
        setColor(new Color(180, 25, 137));

        renderers = new HashMap<>();
        renderers.put(GraphType.ARC, new NestedArcRenderer(this));
        renderers.put(GraphType.PROPORTIONAL_ARC, new ProportionalArcRenderer(this));
        renderers.put(GraphType.BLOCK, new PEBlockRenderer(this));
    }

    private void init(List<BedPEFeature> featureList, Genome genome) {

        this.featureMap = new HashMap<>();

        for (BedPEFeature f : featureList) {
            String key = genome == null ? f.chr1 : genome.getCanonicalChrName(f.chr1);
            if (f.chr1.equals(f.chr2)) {
                addToMap(f, key);
            } else {
                addToMap(new BedPEInterFeature(f, 1), key);
                String key2 = genome == null ? f.chr2 : genome.getCanonicalChrName(f.chr2);
                addToMap(new BedPEInterFeature(f, 2), key2);
            }
        }
        featureMap.put(Globals.CHR_ALL, createWGFeatures(featureList, genome));

        // Sort feature lists by "start" (minimum of start1, start2)
        featureMap.values().forEach(flist -> Collections.sort(flist, (o1, o2) -> o1.getStart() - o2.getStart()));

        // Pack features for block renderer -- TODO do this lazily?
//        for(List<BedPE> flist : featureMap.values()) {
//            BedPEUtils.packFeatures(flist, 100);
//        }
    }

    private void addToMap(BedPE f, String key) {
        List<BedPE> features = featureMap.get(key);
        double maxScore = 0;
        if (features == null) {
            features = new ArrayList<>();
            featureMap.put(key, features);
        }
        features.add(f);

    }


    private List<BedPE> createWGFeatures(List<BedPEFeature> features, Genome genome) {

        List<BedPE> wgFeatures = new ArrayList<>(features.size());

        for (BedPEFeature f : features) {

            int start1 = genome.getGenomeCoordinate(f.chr1, f.start1);
            int end1 = genome.getGenomeCoordinate(f.chr1, f.end1);
            int start2 = genome.getGenomeCoordinate(f.chr2, f.start2);
            int end2 = genome.getGenomeCoordinate(f.chr2, f.end2);
            BedPEFeature wgFeature = new BedPEFeature(Globals.CHR_ALL, start1, end1, Globals.CHR_ALL, start2, end2);

            wgFeature.name = f.name;
            wgFeature.score = f.score;
            wgFeature.thickness = f.thickness;
            wgFeature.color = f.color;
            wgFeature.attributes = f.attributes;


            wgFeatures.add(wgFeature);

        }
        return wgFeatures;
    }

    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {
        return true;
    }

    @Override
    public void load(ReferenceFrame frame) {
        // Nothing to do, this track is pre-loaded
    }

    private List<BedPE> getFeaturesOverlapping(String chr, double start, double end) {

        List<BedPE> features = new ArrayList<>();
        List<BedPE> allFeatures = featureMap.get(chr);

        // TODO - optimize
        if (allFeatures != null) {
            for (BedPE f : allFeatures) {
                if (f.getStart() > end) break;

                if (f.getEnd() >= start) {
                    features.add(f);
                }
            }
        }
        return features;
    }

    @Override
    public void render(RenderContext context, Rectangle trackRectangle) {

        Graphics2D g2d = context.getGraphics();
        Rectangle clip = new Rectangle(g2d.getClip().getBounds());
        g2d.setClip(trackRectangle.intersection(clip.getBounds()));
        context.clearGraphicsCache();


        try {
            String chr = context.getReferenceFrame().getChrName();
            List<BedPE> features = getFeaturesOverlapping(chr, context.getOrigin(), context.getEndLocation());
            if (features != null && features.size() > 0) {
                renderers.get(graphType).render(features, context, trackRectangle);
            }

        } finally {
            context.clearGraphicsCache();
            g2d.setClip(clip);
        }
    }


    @Override
    public IGVPopupMenu getPopupMenu(TrackClickEvent te) {

        IGVPopupMenu menu = new IGVPopupMenu();

        menu.add(TrackMenuUtils.getTrackRenameItem(Collections.singleton(BedPETrack.this)));

        JMenuItem item = new JMenuItem("Set Track Height...");
        item.addActionListener(evt -> TrackMenuUtils.changeTrackHeight(Collections.singleton(BedPETrack.this)));
        menu.add(item);

        item = new JMenuItem("Set Track Color...");
        item.addActionListener(evt -> TrackMenuUtils.changeTrackColor(Collections.singleton(BedPETrack.this)));
        menu.add(item);


        menu.addSeparator();
        menu.add(new JLabel("<html><b>Graph Type</b>"));
        //enum GraphType {BLOCK, ARC, PROPORTIONAL_ARC}
        ButtonGroup group = new ButtonGroup();
        Map<String, GraphType> modes = new LinkedHashMap<>(4);
        modes.put("Nested Arcs", GraphType.ARC);
        modes.put("Proportional Arcs", GraphType.PROPORTIONAL_ARC);
        //modes.put("Blocks", GraphType.BLOCK);

        for (final Map.Entry<String, GraphType> entry : modes.entrySet()) {
            JRadioButtonMenuItem mm = new JRadioButtonMenuItem(entry.getKey());
            mm.setSelected(BedPETrack.this.graphType == entry.getValue());
            mm.addActionListener(evt -> {
                setGraphType(entry.getValue());
                refresh();
            });
            group.add(mm);
            menu.add(mm);
        }

        menu.addSeparator();
        item = new JMenuItem("Toggle Arc Orientation");
        item.addActionListener(evt -> {
            if (direction == UP) {
                direction = Direction.DOWN;
            } else {
                direction = UP;
            }
            IGV.getInstance().repaint();
        });
        menu.add(item);

        item = new JMenuItem("Set Line Thickness...");
        item.addActionListener(e -> {
            String t = MessageUtils.showInputDialog("Line thickness", String.valueOf(thickness));
            if (t != null) {
                try {
                    thickness = Integer.parseInt(t);
                    IGV.getInstance().repaint();
                } catch (NumberFormatException e1) {
                    MessageUtils.showErrorMessage("Line thickness must be an integer", e1);
                }
            }
        });
        menu.add(item);

//        final JCheckBoxMenuItem cbItem = new JCheckBoxMenuItem("Hide Large Features");
//        cbItem.setSelected(hideLargeFeatures);
//        cbItem.addActionListener(e -> {
//            InteractionTrack.this.hideLargeFeatures = cbItem.isSelected();
//            IGV.getInstance().repaint();
//        });
//

        return menu;
    }

    private void setGraphType(GraphType value) {
        // TODO adjust height
        this.graphType = value;
    }

    @Override
    public String getValueStringAt(String chr, double position, int mouseX, int mouseY, ReferenceFrame frame) {
        return super.getValueStringAt(chr, position, mouseX, mouseY, frame);
    }

    @Override
    public void marshalXML(Document document, Element element) {

        super.marshalXML(document, element);

        element.setAttribute("direction", String.valueOf(direction));
        element.setAttribute("thickness", String.valueOf(thickness));
        element.setAttribute("graphType", String.valueOf(graphType));

    }

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);

        if (element.hasAttribute("direction"))
            this.direction = Direction.valueOf(element.getAttribute("direction"));
        if (element.hasAttribute("thickness"))
            this.thickness = Integer.parseInt(element.getAttribute("thickness"));
        if (element.hasAttribute("graphType"))
            this.graphType = GraphType.valueOf(element.getAttribute("graphType").toUpperCase());

    }

}
