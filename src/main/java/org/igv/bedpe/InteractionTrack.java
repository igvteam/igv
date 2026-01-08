package org.igv.bedpe;

import org.igv.Globals;
import org.igv.event.IGVEvent;
import org.igv.event.IGVEventObserver;
import org.igv.hic.HicFile;
import org.igv.jbrowse.CircularViewUtilities;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.renderer.ContinuousColorScale;
import org.igv.renderer.GraphicUtils;
import org.igv.track.AbstractTrack;
import org.igv.track.RenderContext;
import org.igv.track.TrackClickEvent;
import org.igv.track.TrackMenuUtils;
import org.igv.ui.FontManager;
import org.igv.ui.IGV;
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

    private InteractionSource featureSource;
    private JCheckBoxMenuItem autoscaleCB;
    private JMenuItem maxScoreItem;

    InteractionTrack.Direction direction = UP; //DOWN;
    GraphType graphType;  // GraphType.block; //
    private ArcOption arcOption = ArcOption.ALL;
    int thickness = 1;
    boolean autoscale = true;
    double maxScore = -1;
    int gap = 5;
    boolean showBlocks = false;
    boolean useScore = false;
    private Map<GraphType, BedPERenderer> renderers;

    private boolean isHIC;
    ContactMapView contactMapView;
    float transparency = 1.0f;
    String normalization = "NONE";
    private int maxFeatureCount = 5000;

    int[] markerBounds = null;

    transient Map<ReferenceFrame, List<BedPE>> lastRenderedFeatures = new HashMap<>();

    transient Map<ReferenceFrame, LoadedInterval> loadedIntervalMap = new HashMap<>();

    public InteractionTrack() {
    }

    public InteractionTrack(ResourceLocator locator, InteractionSource src) {

        super(locator);

        this.featureSource = src;

        setHeight(250, true);
        setColor(new Color(180, 25, 137));

        renderers = new HashMap<>();
        renderers.put(GraphType.NESTED_ARC, new NestedArcRenderer(this));
        renderers.put(GraphType.PROPORTIONAL_ARC, new ProportionalArcRenderer(this));
        renderers.put(GraphType.BLOCK, new PEBlockRenderer(this));

        this.isHIC = "hic".equals(locator.getFormat());
        if (isHIC) {
            graphType = GraphType.NESTED_ARC;
            useScore = true;
            transparency = 0.1f;
            this.setColor(Color.red);
        } else {

            String typeString = PreferencesManager.getPreferences().get(Constants.ARC_TYPE);
            if (typeString != null) {
                try {
                    graphType = GraphType.valueOf(typeString);
                } catch (IllegalArgumentException e) {
                    log.error("Illegal graph type: " + typeString, e);
                    graphType = GraphType.NESTED_ARC; // default
                }
            } else {
                graphType = GraphType.PROPORTIONAL_ARC;
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
    public void render(RenderContext context, Rectangle trackRectangle) {


        Graphics2D g2d = context.getGraphics();
        Rectangle clip = new Rectangle(g2d.getClip().getBounds());
        g2d.setClip(trackRectangle.intersection(clip.getBounds()));
        context.clearGraphicsCache();

        final ReferenceFrame referenceFrame = context.getReferenceFrame();
        if (!isShowFeatures(referenceFrame)) {
            String message = "Zoom in to see features, or right-click to increase Feature Visibility Window.";
            GraphicUtils.drawCenteredText(message, trackRectangle, context.getGraphics());
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
            List<BedPE> filteredFeatures;
            if (isHIC && interval.zoom < referenceFrame.getZoom()) {
                // In hic mode we limit interactions to those in view plus a margin of one screen width to either side.
                // If zooming in this means we have to filter the features from the previous zoom level that are outside
                // of this range.  Not doing so leads to inconsistent rendering when loading for the current zoom
                // completes and repaints.
                int start = (int) referenceFrame.getOrigin();
                int end = (int) referenceFrame.getEnd();
                int w = (end - start);
                int finalStart = start - w;
                int finalEnd = end + w;
                filteredFeatures = features.stream()
                        .takeWhile(f -> f.getStart() <= finalEnd)
                        .filter(f -> f.getEnd() >= finalStart)
                        .toList();
            } else {
                filteredFeatures = features;
            }


            if (filteredFeatures != null && !filteredFeatures.isEmpty()) {

                if (graphType == GraphType.PROPORTIONAL_ARC) {
                    if (autoscale || maxScore <= 0) {
                        maxScore = autoscale(features);
                    }
                    drawScale(context, trackRectangle);
                }

                renderers.get(graphType).render(filteredFeatures, context, trackRectangle, this.arcOption);
            }
            if (showBlocks) {
                renderers.get(GraphType.BLOCK).render(filteredFeatures, context, trackRectangle, this.arcOption);
            }
            if(contactMapView != null && !FrameManager.isGeneListMode()) {
                //contactMapView.setReferenceFrame(referenceFrame);
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

        if (!this.isHIC) {
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


        final JMenuItem transparencyItem = new JMenuItem("Set Transparency...");
        transparencyItem.addActionListener(e -> {
            final JSlider slider = new JSlider(1, 100, (int) (InteractionTrack.this.transparency * 100));
            slider.setMajorTickSpacing(10);
            slider.setPaintTicks(true);

            // Create a label to show the current value
            final JLabel valueLabel = new JLabel(String.format("%.2f", InteractionTrack.this.transparency));

            slider.addChangeListener(changeEvent -> {
                JSlider source = (JSlider) changeEvent.getSource();
                float value = source.getValue() / 100.0f;
                InteractionTrack.this.transparency = value;
                valueLabel.setText(String.format("%.2f", value));
                InteractionTrack.this.repaint();
            });

            JPanel panel = new JPanel(new BorderLayout());
            panel.add(slider, BorderLayout.CENTER);
            panel.add(valueLabel, BorderLayout.SOUTH);

            final Frame parent = IGV.hasInstance() ? IGV.getInstance().getMainFrame() : null;
            JOptionPane.showMessageDialog(parent, panel, "Set Transparency for " + InteractionTrack.this.getDisplayName(), JOptionPane.PLAIN_MESSAGE);
        });
        menu.add(transparencyItem);

        if (isHIC) {
            final JMenuItem maxFeatureCountItem = new JMenuItem("Set Max Feature Count...");
            maxFeatureCountItem.addActionListener(e -> {
                final JSlider slider = new JSlider(1000, 20000, InteractionTrack.this.maxFeatureCount);
                slider.setMajorTickSpacing(5000);
                slider.setPaintTicks(true);

                final JLabel valueLabel = new JLabel(String.valueOf(InteractionTrack.this.maxFeatureCount));

                slider.addChangeListener(changeEvent -> {
                    JSlider source = (JSlider) changeEvent.getSource();
                    int value = source.getValue();
                    InteractionTrack.this.maxFeatureCount = value;
                    valueLabel.setText(String.valueOf(value));
                    InteractionTrack.this.loadedIntervalMap.clear();
                    InteractionTrack.this.repaint();
                });

                JPanel panel = new JPanel(new BorderLayout());
                panel.add(slider, BorderLayout.CENTER);
                panel.add(valueLabel, BorderLayout.SOUTH);

                final Frame parent = IGV.hasInstance() ? IGV.getInstance().getMainFrame() : null;
                JOptionPane.showMessageDialog(parent, panel, "Set Max Feature Count for " + InteractionTrack.this.getDisplayName(), JOptionPane.PLAIN_MESSAGE);
            });
            menu.add(maxFeatureCountItem);
        }

        // Add normalization options for HiC tracks
        if (isHIC) {
            List<String> normalizationTypes = featureSource.getNormalizationTypes();
            if (normalizationTypes != null && normalizationTypes.size() > 1) {
                menu.addSeparator();
                menu.add(new JLabel("<html><b>Normalization</b>"));
                ButtonGroup normGroup = new ButtonGroup();
                for (String type : normalizationTypes) {
                    String label = normalizationLabels.getOrDefault(type, type);
                    JRadioButtonMenuItem normItem = new JRadioButtonMenuItem(label);
                    normItem.setSelected(type.equals(normalization));
                    normItem.addActionListener(e -> {
                        this.normalization = type;
                        if(contactMapView != null) {
                            contactMapView.setNormalization(type);
                        }
                        InteractionTrack.this.repaint();

                    });
                    normGroup.add(normItem);
                    menu.add(normItem);
                }
            }

            menu.addSeparator();
            JMenuItem mapItem = new JMenuItem("Open Contact Map View");
            mapItem.setEnabled(contactMapView == null && !FrameManager.isGeneListMode());
            mapItem.addActionListener(e -> {
                ReferenceFrame frame = te.getFrame() != null ? te.getFrame() : FrameManager.getDefaultFrame();
                if (contactMapView == null) {
                    ContinuousColorScale colorScale = this.getColorScale();

                    HicFile hicFile = ((HicSource) featureSource).getHicFile();
                    ContactMapView.showPopup(this, hicFile, normalization, frame, colorScale.getMaxColor());
                }
            });
            menu.add(mapItem);

        } else {
            menu.addSeparator();
            menu.add(TrackMenuUtils.getChangeFeatureWindow(Collections.singletonList(this)));
        }

        // Experimental JBrowse.
        if (PreferencesManager.getPreferences().getAsBoolean(Constants.CIRC_VIEW_ENABLED) &&
                CircularViewUtilities.ping() &&
                !isHIC) {
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

        return menu;
    }

    private void setGraphType(GraphType value) {
        // TODO adjust height
        this.graphType = value;
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
        if (this.isHIC) {
            String nviString = ((HicSource) featureSource).getNVIString();
            if (nviString != null) {
                element.setAttribute("nvi", nviString);
            }
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
        if (element.hasAttribute("nvi")) {
            String nviString = element.getAttribute("nvi");
            ((HicSource) featureSource).setNVIString(nviString);
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
