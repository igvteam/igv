package org.broad.igv.bedpe;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.track.AbstractTrack;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.TrackClickEvent;
import org.broad.igv.track.TrackMenuUtils;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.FeatureCache;
import org.broad.igv.util.ResourceLocator;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.*;
import java.util.List;

import static org.broad.igv.bedpe.InteractionTrack.Direction.UP;
import static org.broad.igv.track.TrackMenuUtils.refresh;


/**
 * Created by jrobinso on 6/29/18.
 */
public class InteractionTrack extends AbstractTrack {

    private static Logger log = Logger.getLogger(InteractionTrack.class);

    protected static final int AXIS_AREA_WIDTH = 60;
    protected static Color axisLineColor = new Color(255, 180, 180);
    private JCheckBoxMenuItem autoscaleCB;
    private JMenuItem maxScoreItem;
    private List<BedPE> wgFeatures;


    enum Direction {UP, DOWN}

    enum GraphType {BLOCK, NESTED_ARC, PROPORTIONAL_ARC}

    private Genome genome;
    InteractionTrack.Direction direction = UP; //DOWN;
    GraphType graphType;  // GraphType.block; //
    int thickness = 1;
    boolean autoscale = true;
    double maxScore = -1;
    int gap = 5;
    boolean showBlocks = false;


    //private Map<String, List<BedPE>> featureMap;
    private Map<GraphType, BedPERenderer> renderers;
    private FeatureCache<BedPE> featureCache;

    public InteractionTrack() {
    }

    public InteractionTrack(ResourceLocator locator, BedPEParser.Dataset dataset, Genome genome) {

        super(locator);
        init(dataset.features, genome);
        this.genome = genome;
        setHeight(250, true);
        setColor(new Color(180, 25, 137));

        renderers = new HashMap<>();
        renderers.put(GraphType.NESTED_ARC, new NestedArcRenderer(this));
        renderers.put(GraphType.PROPORTIONAL_ARC, new ProportionalArcRenderer(this));
        renderers.put(GraphType.BLOCK, new PEBlockRenderer(this));

        String typeString = PreferencesManager.getPreferences().get(Constants.ARC_TYPE);
        if (typeString != null) {
            try {
                graphType = GraphType.valueOf(typeString);
            } catch (IllegalArgumentException e) {
                log.error("Illegal graph type: " + typeString, e);
                graphType = GraphType.NESTED_ARC; // default
            }
        } else {
            graphType = dataset.type == BedPEParser.DatasetType.TENX ? GraphType.PROPORTIONAL_ARC : GraphType.NESTED_ARC;
        }


        String directionString = PreferencesManager.getPreferences().get(Constants.ARC_DIRECTION);
        if (directionString != null) {
            try {
                direction = Direction.valueOf(directionString);
            } catch (IllegalArgumentException e) {
                log.error("Illegal arc direction: " + directionString, e);
                direction = UP; // default
            }
        } else {
            direction = UP;
        }

        String blockString = PreferencesManager.getPreferences().get(Constants.ARC_BLOCKS);
        if (blockString != null) {
            try {
                showBlocks = Boolean.valueOf(blockString);
            } catch (IllegalArgumentException e) {
                log.error("Illegal arc blocks option: " + blockString, e);
            }
        }
    }

    private void init(List<BedPEFeature> featureList, Genome genome) {

        List<BedPE> newList = new ArrayList<>((int) (1.2 * featureList.size()));

        for (BedPEFeature f : featureList) {
            String key = genome == null ? f.chr1 : genome.getCanonicalChrName(f.chr1);
            if (f.chr1.equals(f.chr2)) {
                newList.add(f);
            } else {
                newList.add(new BedPEInterFeature(f, 1));
                newList.add(new BedPEInterFeature(f, 2));
            }
        }

        featureCache = new FeatureCache<>(newList, 50);

        wgFeatures = createWGFeatures(featureList, genome);


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

        if(chr.equals(Globals.CHR_ALL)) {
            return wgFeatures;
        } else {
            return featureCache.getFeatures(chr, (int) start, (int) end);
        }
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

                if (graphType == GraphType.PROPORTIONAL_ARC) {
                    if(autoscale || maxScore <= 0) {
                        maxScore = autoscale(features);
                    }
                    drawScale(context, trackRectangle);
                }

                renderers.get(graphType).render(features, context, trackRectangle);
            }
            if (showBlocks) {
                renderers.get(GraphType.BLOCK).render(features, context, trackRectangle);
            }

        } finally {
            context.clearGraphicsCache();
            g2d.setClip(clip);
        }

    }

    /**
     * Draw scale in top left of rectangle

     * @param context
     * @param arect
     */
    public  void drawScale(RenderContext context, Rectangle arect){
        if (context.multiframe == false) {
            Graphics2D g = context.getGraphic2DForColor(Color.black);
            Font font = g.getFont();
            Font smallFont = FontManager.getFont(8);
            try {
                g.setFont(smallFont);
                String minString = "0";
                String fmtString = maxScore > 10 ? "%.0f" : "%.2f";
                String maxString = String.format(fmtString, maxScore);
                String scale = "[" + minString + " - " + maxString + "]";
                g.drawString(scale, arect.x + 5, arect.y + 10);

            } finally {
                g.setFont(font);
            }
        }
    }


    /**
     * Autoscale max height -- specific to proportional arc mode
     *
     * @param features
     */
    private double autoscale(List<BedPE> features) {
        double maxScore = 0;
        for (BedPE f : features) {
            maxScore = Math.max(maxScore, f.getScore());
        }
        return maxScore;
    }


    @Override
    public IGVPopupMenu getPopupMenu(TrackClickEvent te) {

        IGVPopupMenu menu = new IGVPopupMenu();

        menu.add(TrackMenuUtils.getTrackRenameItem(Collections.singleton(InteractionTrack.this)));

        JMenuItem item = new JMenuItem("Set Track Height...");
        item.addActionListener(evt -> TrackMenuUtils.changeTrackHeight(Collections.singleton(InteractionTrack.this)));
        menu.add(item);

        item = new JMenuItem("Set Track Color...");
        item.addActionListener(evt -> TrackMenuUtils.changeTrackColor(Collections.singleton(InteractionTrack.this)));
        menu.add(item);


        menu.addSeparator();
        menu.add(new JLabel("<html><b>Graph Type</b>"));
        //enum GraphType {BLOCK, ARC, PROPORTIONAL_ARC}
        ButtonGroup group = new ButtonGroup();
        Map<String, GraphType> modes = new LinkedHashMap<>(4);
        modes.put("Nested Arcs", GraphType.NESTED_ARC);
        modes.put("Proportional Arcs", GraphType.PROPORTIONAL_ARC);
        //modes.put("Blocks", GraphType.BLOCK);

        for (final Map.Entry<String, GraphType> entry : modes.entrySet()) {
            JRadioButtonMenuItem mm = new JRadioButtonMenuItem(entry.getKey());
            mm.setSelected(InteractionTrack.this.graphType == entry.getValue());
            mm.addActionListener(evt -> {
                setGraphType(entry.getValue());
                PreferencesManager.getPreferences().put(Constants.ARC_TYPE, entry.getValue().toString());
                autoscaleCB.setEnabled(graphType == GraphType.PROPORTIONAL_ARC);
                maxScoreItem.setEnabled(graphType == GraphType.PROPORTIONAL_ARC);
                refresh();
            });
            group.add(mm);
            menu.add(mm);
        }

        menu.addSeparator();
        JCheckBoxMenuItem showBlocksCB = new JCheckBoxMenuItem("Show Blocks");
        showBlocksCB.setSelected(showBlocks);
        showBlocksCB.addActionListener(e -> {
            showBlocks = showBlocksCB.isSelected();
            PreferencesManager.getPreferences().put(Constants.ARC_BLOCKS, String.valueOf(showBlocksCB.isSelected()));
            refresh();
        });
        menu.add(showBlocksCB);

        menu.addSeparator();
        autoscaleCB = new JCheckBoxMenuItem("Autoscale");
        autoscaleCB.setSelected(autoscale);
        autoscaleCB.addActionListener(e -> {
            autoscale = autoscaleCB.isSelected();
            refresh();
        });
        menu.add(autoscaleCB);

        maxScoreItem = new JMenuItem("Set Max Score...");
        maxScoreItem.addActionListener(e -> {
            String maxScoreString = MessageUtils.showInputDialog("Enter maximum score:", String.valueOf(InteractionTrack.this.maxScore));
            if (maxScoreString != null) {
                try {
                    double ms = Double.parseDouble(maxScoreString);
                    if (ms > 0) {
                        maxScore = ms;
                        autoscale = false;
                        refresh();
                    } else {
                        MessageUtils.showMessage("maximum score must be > 0");
                    }
                } catch (NumberFormatException e1) {
                    MessageUtils.showMessage("maximum score must be a number");
                }
            }
        });
        menu.add(maxScoreItem);

        autoscaleCB.setEnabled(graphType == GraphType.PROPORTIONAL_ARC);
        maxScoreItem.setEnabled(graphType == GraphType.PROPORTIONAL_ARC);


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

        List<BedPE> candidates = getFeaturesOverlapping(frame.getChrName(), (int) position, (int) position + 1);

        // Sort candidate features smallest to largest
        Comparator<BedPE> sorter = graphType == GraphType.PROPORTIONAL_ARC ?
                (o1, o2) -> {
                    double score1 = o1.getScore();
                    double score2 = o2.getScore();
                    if (score1 > score2) return 1;
                    else if (score1 < score2) return -1;
                    else return 0;
                } :
                (o1, o2) -> {
                    double d1 = o1.getCenterDistance();
                    double d2 = o2.getCenterDistance();
                    if (d1 > d2) return 1;
                    else if (d1 < d2) return -1;
                    else return 0;
                };

        Collections.sort(candidates, sorter);

        for (BedPE f : candidates) {
            BedPEShape s = f.getShape();
            if (s != null && s.contains(mouseX, mouseY)) {
                return f.getValueString();
            }
        }

        return super.getValueStringAt(chr, position, mouseX, mouseY, frame);
    }

    @Override
    public void marshalXML(Document document, Element element) {

        super.marshalXML(document, element);

        element.setAttribute("direction", String.valueOf(direction));
        element.setAttribute("thickness", String.valueOf(thickness));
        element.setAttribute("graphType", String.valueOf(graphType));
        element.setAttribute("showBlocks", String.valueOf(showBlocks));
        element.setAttribute("autoscale", String.valueOf(autoscale));
        if (!autoscale) {
            element.setAttribute("maxScore", String.valueOf(maxScore));
        }

    }

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);

        if (element.hasAttribute("direction"))
            this.direction = Direction.valueOf(element.getAttribute("direction"));
        if (element.hasAttribute("thickness"))
            this.thickness = Integer.parseInt(element.getAttribute("thickness"));
        if (element.hasAttribute("graphType")) {
            String typeString = element.getAttribute("graphType").toUpperCase();
            if (typeString.equals("ARC")) typeString = "NESTED_ARC";  // backward compatibility
            this.graphType = GraphType.valueOf(typeString);
        }
        if (element.hasAttribute("showBlocks"))
            this.showBlocks = Boolean.parseBoolean(element.getAttribute("showBlocks"));
        if (element.hasAttribute("autoscale")) {
            this.autoscale = Boolean.parseBoolean(element.getAttribute("autoscale"));
        }
        if (element.hasAttribute("maxScore")) {
            this.maxScore = Double.parseDouble(element.getAttribute("maxScore"));
        }

    }

}
