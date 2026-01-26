package org.igv.bedpe;

import org.igv.Globals;
import org.igv.event.IGVEvent;
import org.igv.event.IGVEventObserver;
import org.igv.jbrowse.CircularViewUtilities;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.renderer.GraphicUtils;
import org.igv.track.AbstractTrack;
import org.igv.track.RenderContext;
import org.igv.track.TrackClickEvent;
import org.igv.track.TrackMenuUtils;
import org.igv.ui.FontManager;
import org.igv.ui.panel.FrameManager;
import org.igv.ui.panel.IGVPopupMenu;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.ui.util.MessageUtils;
import org.igv.ui.util.UIUtilities;
import org.igv.util.ResourceLocator;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import javax.swing.*;
import java.awt.*;
import java.io.IOException;
import java.util.*;
import java.util.List;

import static org.igv.bedpe.InteractionTrack.Direction.UP;

/**
 * Created by jrobinso on 6/29/18.
 */
public class InteractionTrack extends AbstractTrack implements IGVEventObserver {

    private static final Logger log = LogManager.getLogger(InteractionTrack.class);

    enum Direction {UP, DOWN}

    enum GraphType {BLOCK, NESTED_ARC, PROPORTIONAL_ARC}

    enum ArcOption {ALL, ONE_END, BOTH_ENDS}

    static Map<String, String> normalizationLabels = new LinkedHashMap<>();

    static {
        normalizationLabels.put("NONE", "None");
        normalizationLabels.put("VC", "Coverage");
        normalizationLabels.put("VC_SQRT", "Coverage - Sqrt");
        normalizationLabels.put("KR", "Balanced (Knight-Ruiz)");
        normalizationLabels.put("INTER_VC", "Interchromosomal Coverage");
        normalizationLabels.put("INTER_VC_SQRT", "Interchromosomal Coverage - Sqrt");
        normalizationLabels.put("INTER_KR", "Interchromosomal Balanced");
        normalizationLabels.put("GW_VC", "Genome-wide Coverage");
        normalizationLabels.put("GW_VC_SQRT", "Genome-wide Coverage - Sqrt");
        normalizationLabels.put("GW_KR", "Genome-wide Balanced");
    }

    protected InteractionSource featureSource;
    private JCheckBoxMenuItem autoscaleCB;
    private JMenuItem maxScoreItem;

    InteractionTrack.Direction direction = UP; //DOWN;
    protected GraphType graphType = GraphType.NESTED_ARC;  // GraphType.block; //
    private ArcOption arcOption = ArcOption.ALL;
    int thickness = 1;
    boolean autoscale = true;
    double maxScore = -1;
    int gap = 5;
    boolean showBlocks = false;
    protected boolean isHIC = false;
    private Map<GraphType, BedPERenderer> renderers;

    ContactMapView contactMapView;
    float transparency = 1.0f;
    protected String normalization = "NONE";
    protected int maxFeatureCount = 20000;

    int[] markerBounds = null;

    transient Map<ReferenceFrame, List<BedPE>> lastRenderedFeatures = new HashMap<>();

    transient Map<ReferenceFrame, LoadedInterval> loadedIntervalMap = new HashMap<>();

    public InteractionTrack() {
    }

    public InteractionTrack(ResourceLocator locator, InteractionSource src) {

        super(locator);

        this.featureSource = src;

        setHeight(250);
        setDefaultColor( new Color(180, 25, 137));

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
            }
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
                showBlocks = Boolean.parseBoolean(blockString);
            } catch (IllegalArgumentException e) {
                log.error("Illegal arc blocks option: " + blockString, e);
            }
        }
    }

    void setContactMapView(ContactMapView contactMapView) {
        this.contactMapView = contactMapView;
    }

    protected boolean isShowFeatures(ReferenceFrame frame) {

        if (frame.getChrName().equals(Globals.CHR_ALL)) {
            return true;
        } else {
            double windowSize = frame.getEnd() - frame.getOrigin();
            int vw = getVisibilityWindow();
            return (vw <= 0 || windowSize <= vw);
        }
    }

    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {

        LoadedInterval interval = loadedIntervalMap.get(frame);
        final boolean b = interval != null && interval.contains(frame.getChrName(),
                (int) frame.getOrigin(),
                (int) frame.getEnd(),
                frame.getZoom(),
                normalization);
        return b;
    }

    @Override
    public void load(ReferenceFrame frame) {
        String chr = frame.getChrName();
        int start = (int) frame.getOrigin();
        int end = (int) frame.getEnd();
        int zoom = frame.getZoom();

        // Expand region by half width to enable panning without reloading
        int w = (end - start) / 2;
        start = Math.max(0, start - w);
        end = end + w;

        try {
            List<BedPE> features = featureSource.getFeatures(chr, start, end, frame.getScale(), normalization, maxFeatureCount);
            LoadedInterval interval = new LoadedInterval(chr, start, end, zoom, normalization, features);
            loadedIntervalMap.put(frame, interval);
        } catch (IOException e) {
            log.error("Error loading features", e);
        }
    }

    @Override
    public void render(RenderContext context, Rectangle visibleRect) {


        Graphics2D g2d = context.getGraphics();
        Rectangle clip = new Rectangle(g2d.getClip().getBounds());
        g2d.setClip(visibleRect.intersection(clip.getBounds()));
        context.clearGraphicsCache();

        final ReferenceFrame referenceFrame = context.getReferenceFrame();
        if (!isShowFeatures(referenceFrame)) {
            String message = "Zoom in to see features, or right-click to increase Feature Visibility Window.";
            GraphicUtils.drawCenteredText(message, visibleRect, context.getGraphics());
            return;
        }

        try {
            String chr = referenceFrame.getChrName();

            if (normalization != null && !normalization.equals("NONE") &&
                    !featureSource.hasNormalizationVector(normalization, chr, context.getScale())) {
                String message = "Normalization '" + normalization + "' not available at this resolution. Switching normalization to 'NONE'.";
                normalization = "NONE";
                UIUtilities.invokeOnEventThread(() -> MessageUtils.showMessage(message));
            }

            final LoadedInterval interval = loadedIntervalMap.get(referenceFrame);
            if (interval == null) {
                return;
            }

            List<BedPE> features = interval.features();
            List<BedPE> filteredFeatures = filterFeaturesForZoom(features, interval, referenceFrame);


            if (filteredFeatures != null && !filteredFeatures.isEmpty()) {

                if (graphType == GraphType.PROPORTIONAL_ARC) {
                    if (autoscale || maxScore <= 0) {
                        maxScore = autoscale(features);
                    }
                    drawScale(context, visibleRect);
                }

                renderers.get(graphType).render(filteredFeatures, context, visibleRect, this.arcOption);
            }
            if (showBlocks) {
                renderers.get(GraphType.BLOCK).render(filteredFeatures, context, visibleRect, this.arcOption);
            }
            if (contactMapView != null && !FrameManager.isGeneListMode()) {
            }

        } finally {
            context.clearGraphicsCache();
            g2d.setClip(clip);
        }
    }

    /**
     * Draw scale in top left of rectangle
     *
     * @param context
     * @param arect
     */
    public void drawScale(RenderContext context, Rectangle arect) {
        if (!context.multiframe) {
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

    /**
     * Hook method for subclasses to filter features based on zoom level.
     * Default implementation returns features unchanged.
     */
    protected List<BedPE> filterFeaturesForZoom(List<BedPE> features, LoadedInterval interval, ReferenceFrame referenceFrame) {
        return features;
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

        item = new JMenuItem("Unset Track Color");
        item.addActionListener(evt -> {
            this.setColor(null);
            repaint();
        });
        menu.add(item);

        if (!isHIC) {
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
                    this.graphType = entry.getValue();
                    PreferencesManager.getPreferences().put(Constants.ARC_TYPE, entry.getValue().toString());
                    autoscaleCB.setEnabled(graphType == GraphType.PROPORTIONAL_ARC);
                    maxScoreItem.setEnabled(graphType == GraphType.PROPORTIONAL_ARC);
                    repaint();
                });
                group.add(mm);
                menu.add(mm);
            }

            menu.addSeparator();
            menu.add(new JLabel("<html><b>Arcs</b>"));
            ButtonGroup group2 = new ButtonGroup();
            Map<String, ArcOption> modes2 = new LinkedHashMap<>(4);
            modes2.put("All", ArcOption.ALL);
            modes2.put("One End In View", ArcOption.ONE_END);
            modes2.put("Both Ends In View", ArcOption.BOTH_ENDS);
            //modes.put("Blocks", GraphType.BLOCK);

            for (final Map.Entry<String, ArcOption> entry : modes2.entrySet()) {
                JRadioButtonMenuItem mm = new JRadioButtonMenuItem(entry.getKey());
                mm.setSelected(InteractionTrack.this.arcOption == entry.getValue());
                mm.addActionListener(evt -> {
                    InteractionTrack.this.arcOption = (entry.getValue());
                    repaint();
                });
                group2.add(mm);
                menu.add(mm);
            }
        }

        menu.addSeparator();
        JCheckBoxMenuItem showBlocksCB = new JCheckBoxMenuItem("Show Blocks");
        showBlocksCB.setSelected(showBlocks);
        showBlocksCB.addActionListener(e -> {
            showBlocks = showBlocksCB.isSelected();
            PreferencesManager.getPreferences().put(Constants.ARC_BLOCKS, String.valueOf(showBlocksCB.isSelected()));
            repaint();
        });
        menu.add(showBlocksCB);

        if (!isHIC) {
            menu.addSeparator();
            autoscaleCB = new JCheckBoxMenuItem("Autoscale");
            autoscaleCB.setSelected(autoscale);
            autoscaleCB.addActionListener(e -> {
                autoscale = autoscaleCB.isSelected();
                repaint();
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
                            repaint();
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
        }

        menu.addSeparator();
        item = new JMenuItem("Toggle Arc Orientation");
        item.addActionListener(evt -> {
            if (direction == UP) {
                direction = Direction.DOWN;
            } else {
                direction = UP;
            }
            repaint();
        });
        menu.add(item);

        item = new JMenuItem("Set Line Thickness...");
        item.addActionListener(e -> {
            String t = MessageUtils.showInputDialog("Line thickness", String.valueOf(thickness));
            if (t != null) {
                try {
                    thickness = Integer.parseInt(t);
                    repaint();
                } catch (NumberFormatException e1) {
                    MessageUtils.showErrorMessage("Line thickness must be an integer", e1);
                }
            }
        });
        menu.add(item);


        if (isHIC) {
            addHICItems(te, menu);

        } else {
            // Not hic
            menu.addSeparator();
            menu.add(TrackMenuUtils.getChangeFeatureWindow(Collections.singletonList(this)));


            // Experimental JBrowse.
            if (PreferencesManager.getPreferences().getAsBoolean(Constants.CIRC_VIEW_ENABLED) &&
                    CircularViewUtilities.ping()) {
                menu.addSeparator();
                JMenuItem circViewItem = new JMenuItem("Add Features to Circular View");
                circViewItem.addActionListener(e -> {
                    List<ReferenceFrame> frames = te.getFrame() != null ?
                            Collections.singletonList(te.getFrame()) :
                            FrameManager.getFrames();
                    List<? extends BedPE> visibleFeatures = getVisibleFeatures(frames);
                    CircularViewUtilities.sendBedpeToJBrowse(visibleFeatures, InteractionTrack.this.getName(), InteractionTrack.this.getColor());
                });
                menu.add(circViewItem);
                menu.addSeparator();
            }
        }

        return menu;
    }

    void addHICItems(TrackClickEvent te, IGVPopupMenu menu) {
        // Override in HIC subclass
    }

    @Override
    public void setColor(Color color) {
        super.setColor(color);
    }

    @Override
    public String getValueStringAt(String chr, double position, int mouseX, int mouseY, ReferenceFrame frame) {

        // Expand range a little bit -- this should be done in pixels
        double tolerance = frame.getScale() * 3;

        LoadedInterval interval = loadedIntervalMap.get(frame);
        if (interval == null) return "";

        List<BedPE> features = interval.features();
        if (features == null) return "";

        List<BedPE> candidates = new ArrayList<>();
        for (BedPE bedPE : features) {
            if (bedPE.getEnd() < position - tolerance) continue;
            if (bedPE.getStart() > position + tolerance) break;
            candidates.add(bedPE);
        }

        // Sort candidate features smallest to largest
        Comparator<BedPE> sorter = graphType == GraphType.PROPORTIONAL_ARC ?
                Comparator.comparingDouble(BedPE::getScore) :
                Comparator.comparingDouble(BedPE::getCenterDistance);

        candidates.sort(sorter);

        for (BedPE f : candidates) {
            BedPEShape s = f.getShape();
            if (s == null || s.contains(mouseX, mouseY)) {
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
        element.setAttribute("arcOption", String.valueOf(arcOption));
        element.setAttribute("showBlocks", String.valueOf(showBlocks));
        element.setAttribute("autoscale", String.valueOf(autoscale));
        if (transparency != 1.0f) {
            element.setAttribute("transparency", String.valueOf(transparency));
        }
        if (!autoscale) {
            element.setAttribute("maxScore", String.valueOf(maxScore));
        }

    }

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);

        if (element.hasAttribute("arcOption")) {
            String typeString = element.getAttribute("arcOption").toUpperCase();
            this.arcOption = ArcOption.valueOf(typeString);
        }
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
        if (element.hasAttribute("transparency")) {
            this.transparency = Float.parseFloat(element.getAttribute("transparency"));
        }
    }


    /**
     * Return features visible in the supplied frames
     */
    /**
     * Return features visible in the supplied frames
     */
    public List<BedPE> getVisibleFeatures(List<ReferenceFrame> frames) {

        if (frames.isEmpty()) {
            return Collections.emptyList();
        } else {
            List<BedPE> inView = new ArrayList<>();
            for (ReferenceFrame f : frames) {
                LoadedInterval interval = loadedIntervalMap.get(f);
                if (interval != null) {
                    inView.addAll(interval.features());
                }
            }
            return inView;
        }
    }

    @Override
    public void receiveEvent(IGVEvent event) {
        if (event instanceof FrameManager.ChangeEvent) {
            Set<ReferenceFrame> frames = new HashSet<>(FrameManager.getFrames());

            // Remove cached intervals for frames that are no longer present
            loadedIntervalMap.keySet().removeIf(f -> !frames.contains(f));
        }
    }

    public void setMarkerBounds(int[] markerBounds) {
        this.markerBounds = markerBounds;
        this.repaint();
    }


    public record LoadedInterval(String chr, int start, int end, int zoom, String normalization,
                                 List<BedPE> features) {

        String getKey() {
            return chr + "_" + start + "_" + end + "_" + zoom + "_" + normalization;
        }

        boolean contains(String chr, int start, int end, int zoom, String normalization) {
            return this.chr.equals(chr) &&
                    this.start <= start &&
                    this.end >= end &&
                    (this.zoom == -1 || this.zoom == zoom) &&
                    (this.normalization == null || this.normalization.equals(normalization));
        }
    }
}
