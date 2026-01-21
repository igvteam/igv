package org.igv.track;


import htsjdk.tribble.Feature;
import org.igv.Globals;
import org.igv.event.IGVEventBus;
import org.igv.event.IGVEventObserver;
import org.igv.feature.IGVFeature;
import org.igv.feature.Strand;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.renderer.*;
import org.igv.session.SessionAttribute;
import org.igv.ui.FontManager;
import org.igv.ui.IGV;
import org.igv.ui.TooltipTextFrame;
import org.igv.ui.color.ColorUtilities;
import org.igv.ui.panel.AttributeHeaderPanel;
import org.igv.ui.panel.IGVPopupMenu;
import org.igv.ui.panel.MouseableRegion;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.ui.util.UIUtilities;
import org.igv.util.ResourceLocator;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.*;
import java.util.List;

import static org.igv.prefs.Constants.*;

/**
 * @author jrobinso
 */

public abstract class AbstractTrack implements Track {

    private static Logger log = LogManager.getLogger(AbstractTrack.class);

    public static final DisplayMode DEFAULT_DISPLAY_MODE = DisplayMode.COLLAPSED;
    public static final int DEFAULT_HEIGHT = -1;
    public static final int VISIBILITY_WINDOW = -1;
    public static final boolean DEFAULT_SHOW_FEATURE_NAMES = true;
    protected final boolean darkMode;

    private ResourceLocator resourceLocator;
    private String id;
    private String sampleId;
    private String attributeKey;
    private String name;
    private String featureInfoURL;
    private boolean itemRGB = true;

    private boolean useScore;
    protected boolean autoScale = true;   // By default, can be overriden by track line
    private float viewLimitMin = 0;     // From UCSC track line
    private float viewLimitMax = 1000;  // From UCSC track line
    protected int fontSize;
    private boolean showDataRange = true;

    private int top;
    protected int minimumHeight = -1;
    private TrackType trackType = TrackType.OTHER;
    private boolean selected = false;
    private boolean visible = true;
    private boolean sortable = true;

    boolean overlaid;
    boolean drawYLine = false;
    float yLine = 0;


    private Map<String, String> attributes = new HashMap();

    private ContinuousColorScale colorScale;

    String autoscaleGroup;

    protected Color defaultColor;
    protected Color color = null;
    protected Color altColor = null;

    protected int visibilityWindow = VISIBILITY_WINDOW;
    private DisplayMode displayMode = DEFAULT_DISPLAY_MODE;

    protected Integer height = DEFAULT_HEIGHT;

    protected DataRange dataRange;
    private boolean showFeatureNames = DEFAULT_SHOW_FEATURE_NAMES;
    private String trackLine = null;

    public AbstractTrack() {
        this.darkMode = Globals.isDarkMode();
        defaultColor = darkMode ? Color.CYAN : Color.blue.darker();
    }

    public AbstractTrack(
            ResourceLocator dataResourceLocator,
            String id,
            String name) {
        this();
        this.resourceLocator = dataResourceLocator;
        this.id = id;
        this.name = name;
        this.attributeKey = this.name;
        init();
    }

    public AbstractTrack(ResourceLocator locator) {
        this(locator, locator != null ? locator.getPath() : null, locator != null ? locator.getTrackName() : null);
    }

    private void init() {

        try {
            fontSize = PreferencesManager.getPreferences().getAsInt(DEFAULT_FONT_SIZE);
        } catch (Exception e) {
            log.error("Error initializing font size. ", e);
        }

        showDataRange = PreferencesManager.getPreferences().getAsBoolean(CHART_SHOW_DATA_RANGE);
        if (PreferencesManager.getPreferences().getAsBoolean(EXPAND_FEAUTRE_TRACKS)) {
            displayMode = DisplayMode.EXPANDED;
        }
    }

    public void setRendererClass(Class rc) {
        // Ignore by default
    }

    public String getFeatureInfoURL() {
        return featureInfoURL;
    }

    public void setFeatureInfoURL(String featureInfoURL) {
        this.featureInfoURL = featureInfoURL;
    }

    public void setUseScore(boolean useScore) {
        this.useScore = useScore;
    }

    public String getId() {
        return id;
    }


    public void setName(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    public String getDisplayName() {

        String sampleKey = IGV.getInstance().getSession().getTrackAttributeName();
        if (sampleKey != null && sampleKey.trim().length() > 0) {
            String name = getAttributeValue(sampleKey.trim());
            if (name != null) {
                return name;
            }
        }
        return getName();
    }


    public void setSampleId(String sampleId) {
        this.sampleId = sampleId;
    }

    public void renderName(Graphics2D g2D, Rectangle trackRectangle, Rectangle visibleRectangle) {

        Rectangle rect = getDisplayableRect(trackRectangle, visibleRectangle);

        String trackName = getDisplayName();
        if ((trackName != null)) {

            if (rect.getHeight() > 3) {

                // Calculate fontsize
                int gap = Math.min(4, rect.height / 3);
                int fs = Math.min(fontSize, rect.height - gap);

                Font font = FontManager.getFont(fs);
                g2D.setFont(font);

                GraphicUtils.drawWrappedText(trackName, rect, g2D, false);

                //g2D.dispose();
            }
        }
    }


    public void renderAttributes(Graphics2D graphics, Rectangle trackRectangle, Rectangle visibleRect,
                                 List<String> names, List<MouseableRegion> mouseRegions) {

        int x = trackRectangle.x;

        for (String name : names) {
            String key = name.toUpperCase();
            String attributeValue = getAttributeValue(key);
            if (attributeValue != null) {
                Rectangle rect = new Rectangle(x, trackRectangle.y, AttributeHeaderPanel.ATTRIBUTE_COLUMN_WIDTH,
                        trackRectangle.height);
                graphics.setColor(AttributeManager.getInstance().getColor(key, attributeValue));
                graphics.fill(rect);
                mouseRegions.add(new MouseableRegion(rect, key, attributeValue));

                // Border
//                graphics.setColor(Color.lightGray);
//                graphics.fillRect(x + AttributeHeaderPanel.ATTRIBUTE_COLUMN_WIDTH, trackRectangle.y,
//                        AttributeHeaderPanel.COLUMN_BORDER_WIDTH, trackRectangle.height);

            }
            x += AttributeHeaderPanel.ATTRIBUTE_COLUMN_WIDTH + AttributeHeaderPanel.COLUMN_BORDER_WIDTH;

        }
    }


    private Rectangle getDisplayableRect(Rectangle trackRectangle, Rectangle visibleRect) {
        Rectangle rect = null;
        if (visibleRect != null) {
            Rectangle intersectedRect = trackRectangle.intersection(visibleRect);
            if (intersectedRect.height > 15) {
                rect = intersectedRect;
            } else {
                rect = new Rectangle(trackRectangle);
            }
        }
        return rect;

    }


    /**
     * Called to overlay a track on another, presumably previously rendered,
     * track. The default behavior is to do nothing.
     *
     * @param context
     * @param rect
     */
    public void overlay(RenderContext context, Rectangle rect) {
    }

    public void setColor(Color color) {
        this.color = color;
    }

    public Color getColor() {
        return color;
    }

    @Override
    public void setAltColor(Color color) {
        altColor = color;
    }

    public Color getAltColor() {
        return altColor == null ? getColor() : altColor;
    }

    public void setDefaultColor(Color color) {
        this.defaultColor = color;
    }

    @Override
    public Color getDefaultColor() {
        return defaultColor;
    }

    public ResourceLocator getResourceLocator() {
        return resourceLocator;
    }

    public Collection<ResourceLocator> getResourceLocators() {
        return Arrays.asList(getResourceLocator());
    }

    /**
     * Add an attribute to this track and register the key with the attribute panel.
     * <p/>
     * Note:  Attribute keys are case insensitive.  Currently this is implemented
     * by forcing all keys to upper case
     *
     * @param name
     * @param value
     */
    public void setAttributeValue(String name, String value) {
        String key = name.toUpperCase();
        if (key.equals(AttributeManager.GROUP_AUTOSCALE)) {
            autoscaleGroup = value;
        } else {
            attributes.put(key, value);
        }
        AttributeManager.getInstance().addAttribute(getSample(), name, value);

    }

    public void removeAttribute(String name) {
        String key = name.toUpperCase();
        if (key.equals(AttributeManager.GROUP_AUTOSCALE)) {
            autoscaleGroup = null;
        } else {
            attributes.remove(key);
        }
        AttributeManager.getInstance().removeAttribute(getSample(), name);
    }


    /**
     * Return the attribute value.  Attribute lookup occurs in the following order, if all fail null is returned.
     * <p/>
     * (1) the track attribute table
     * (2) by sampleId, as set in the Resource element of a session or load-from-server menu
     * (3) by attributeKey, set from the original track name
     * (4) by full path to the file associated with this track
     *
     * @param attributeName
     * @return
     */
    public String getAttributeValue(String attributeName) {

        String key = attributeName.toUpperCase();
        String value;
        if (key.equals(AttributeManager.GROUP_AUTOSCALE)) {
            value = autoscaleGroup;
            if (value == null) {
                value = getFromAttributeManager(key);
                autoscaleGroup = value;
            }
        } else {
            value = attributes.get(key);
            if (value == null) {
                value = getFromAttributeManager(key);
            }
        }
        return value;
    }

    private String getFromAttributeManager(String key) {
        final AttributeManager attributeManager = AttributeManager.getInstance();
        String value = null;
        if (value == null && getSample() != null) {
            value = attributeManager.getAttribute(getSample(), key);
        }
        if (value == null) {
            value = attributeManager.getAttribute(this.attributeKey, key);
        }
        if (value == null && getResourceLocator() != null && getResourceLocator().getPath() != null) {
            value = attributeManager.getAttribute(getResourceLocator().getPath(), key);
        }
        return value;
    }

    public String getSample() {
        if (sampleId != null) {
            return sampleId;    // Explicitly set sample ID (e.g. from server load XML)
        }
        sampleId = AttributeManager.getInstance().getSampleFor(getName());
        return sampleId != null ? sampleId : getName();
    }


    /**
     * Returns the default height based on the default renderer for the data
     * type, as opposed to the actual renderer in use.  This is done to prevent
     * the track size from changing if renderer is changed.
     *
     * @return
     */
    private int getDefaultHeight() {
        if (getDefaultRenderer() instanceof XYPlotRenderer) {
            return PreferencesManager.getPreferences().getAsInt(CHART_TRACK_HEIGHT_KEY);
        } else {
            return PreferencesManager.getPreferences().getAsInt(TRACK_HEIGHT_KEY);
        }
    }


    /**
     * Returns the default minimum height based on the actual renderer for this track.  Heatmaps default
     * to 1,  all other renderers to 5.
     *
     * @return
     */
    public int getDefaultMinimumHeight() {
        Renderer r = getRenderer();
        if (r != null && HeatmapRenderer.class.isAssignableFrom(r.getClass())) {
            return 1;
        } else {
            return 10;
        }
    }

    public void setMinimumHeight(int minimumHeight) {
        this.minimumHeight = minimumHeight;
    }

    /**
     * Return the actual minimum height if one has been set, otherwise get the default for the current renderer.
     *
     * @return
     */
    public int getMinimumHeight() {
        return minimumHeight < 0 ? getDefaultMinimumHeight() : minimumHeight;
    }

    @Override
    public void setTrackType(TrackType type) {
        this.trackType = type;
    }

    public TrackType getTrackType() {
        return trackType;
    }

    public boolean isVisible() {

        if (visible && getTrackType() == TrackType.MUTATION) {
            // Special rules for mutations.  If display as overlays == true, only show if not overlaid on another
            // track and "show orphaned" is true
            boolean displayOverlays = IGV.getInstance().getSession().getOverlayMutationTracks();
            if (displayOverlays) {
                if (overlaid) {
                    return false;
                } else {
                    return PreferencesManager.getPreferences().getAsBoolean(SHOW_ORPHANED_MUTATIONS);
                }
            }
        }
        return visible;
    }

    public void setVisible(boolean visible) {
        if (this.visible != visible) {
            this.visible = visible;
            if (IGV.hasInstance()) IGV.getInstance().getMainPanel().revalidate();
        }
    }

    public void setOverlayed(boolean bool) {
        this.overlaid = bool;
    }

    public void setSelected(boolean selected) {
        this.selected = selected;
    }

    public boolean isSelected() {
        return selected;
    }

    public void setHeight(int height) {
        setHeight(height, false);
    }

    @Override
    public void setHeight(int preferredHeight, boolean force) {
        if (height < getHeight()) {
            if ((this.getDisplayMode() == DisplayMode.EXPANDED) && (getTrackType() != TrackType.GENE)) {
                this.setDisplayMode(DisplayMode.SQUISHED);
            }
        }

        if (force) {
            this.height = preferredHeight;
        } else {
            this.height = Math.max(getMinimumHeight(), preferredHeight);
        }
    }

    public int getHeight() {
        return (height < 0) ? getDefaultHeight() : height;
    }

    public DataRange getDataRange() {
        if (dataRange == null) {
            // Use the color scale if there is one
            float min = (float) (colorScale == null ? 0 : colorScale.getMinimum());
            float max = (float) (colorScale == null ? 10 : colorScale.getMaximum());
            float baseline = (float) (colorScale == null ? 0 : (colorScale.getNegStart() + colorScale.getPosStart()) / 2);

            setDataRange(new DataRange(min, baseline, max));
        }
        return dataRange;
    }


    public void setDataRange(DataRange axisDefinition) {
        this.dataRange = axisDefinition;
    }

    protected Renderer getDefaultRenderer() {
        TrackType trackType = getTrackType();
        switch (trackType) {
            case RNAI:
            case COPY_NUMBER:
            case ALLELE_SPECIFIC_COPY_NUMBER:
            case GENE_EXPRESSION:
            case DNA_METHYLATION:
            case LOH:
            case CHIP_CHIP:
                return new HeatmapRenderer();
            default:
                return new BarChartRenderer();
        }
    }

    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return new ArrayList();
    }

    public boolean handleDataClick(TrackClickEvent te) {

        if (IGV.getInstance().isShowDetailsOnClick()) {
            return openTooltipWindow(te);
        }
        return false;
    }

    protected boolean openTooltipWindow(TrackClickEvent e) {
        ReferenceFrame frame = e.getFrame();
        final MouseEvent me = e.getMouseEvent();
        String popupText = getValueStringAt(frame.getChrName(), e.getChromosomePosition(), e.getMouseEvent().getX(), e.getMouseEvent().getY(), frame);

        if (popupText != null) {
            Color color = IGV.getInstance().getRootPane().getJMenuBar().getForeground();
            String htmlColor = String.format("#%02x%02x%02x", color.getRed(), color.getGreen(), color.getBlue());
            popupText = "<div style=\"color: " + htmlColor + "\">" + popupText + "</div>";

            final TooltipTextFrame tf = new TooltipTextFrame(getName(), popupText);
            Point p = me.getComponent().getLocationOnScreen();
            tf.setLocation(Math.max(0, p.x + me.getX() - 150), Math.max(0, p.y + me.getY() - 150));

            UIUtilities.invokeOnEventThread(() -> tf.setVisible(true));
            return true;
        }
        return false;
    }

    public void handleNameClick(MouseEvent e) {
        // Do nothing
    }

    @Override
    public void setAutoScale(boolean autoScale) {
        this.autoScale = autoScale;
    }

    /**
     * Set some properties of this track,  usually from a "track line" specification.
     * <p/>
     * TODO -- keep the properties object, rather than copy all the values.
     *
     * @param properties
     */
    public void setProperties(TrackProperties properties) {

        if (properties.isItemRGB() != null) {
            this.itemRGB = properties.isItemRGB();
        }
        if (properties.isUseScore() != null) {
            this.useScore = properties.isUseScore();
        }
        if (properties.getyLine() != null) {
            this.yLine = properties.getyLine();
        }
        this.drawYLine = properties.isDrawYLine();
        this.sortable = properties.isSortable();

        // The viewLimit properties com from UCSC track lines and control "useScore" shading. They
        // are somewhat redundant with dataRange,  but we keep both for now to avoid breaking existing behavior.
        if (properties.getMinValue() != null) {
            this.viewLimitMin = properties.getMinValue();
        }
        if (properties.getMaxValue() != null) {
            this.viewLimitMax = properties.getMaxValue();
        }

        if (properties.getMaxValue() != null) {

            float max = properties.getMaxValue();

            Float minVal = properties.getMinValue();
            float min = (minVal == null) ? 0f : minVal;

            Float midVal = properties.getMidValue();
            float mid = (midVal == null) ? (min < 0 ? 0f : min) : midVal;

            DataRange dr = new DataRange(min, mid, max);
            setDataRange(dr);

            if (properties.isLogScale()) {
                dr.setType(DataRange.Type.LOG);
            }

            // If the user has explicity set a data range and colors apply to heatmap as well
            Color maxColor = properties.getColor();
            Color minColor = properties.getAltColor();
            if (maxColor != null && minColor != null) {

                Float tmp = properties.getNeutralFromValue();
                float neutralFrom = tmp == null ? mid : tmp;
                tmp = properties.getNeutralToValue();
                float neutralTo = tmp == null ? mid : tmp;

                Color midColor = properties.getMidColor();
                if (midColor == null) {
                    midColor = Color.white;
                }
                colorScale = new ContinuousColorScale(neutralFrom, min, neutralTo, max, minColor, midColor, maxColor);
            }
        }

        if (properties.getDisplayMode() != null) {
            this.setDisplayMode(properties.getDisplayMode());
        }

        if (properties.getName() != null) {
            name = properties.getName();
        }
        if (properties.getColor() != null) {
            setColor(properties.getColor());
        }
        if (properties.getAltColor() != null) {
            setAltColor(properties.getAltColor());
        }
        if (properties.getMidColor() != null) {
            //setMidColor(trackProperties.getMidColor());
        }
        if (properties.getHeight() != null) {
            setHeight(properties.getHeight());
        }
        if (properties.getMinHeight() != null) {
            setMinimumHeight(properties.getMinHeight());
        }
        if (properties.getRendererClass() != null) {
            setRendererClass(properties.getRendererClass());
            if (properties.getRendererClass() == PointsRenderer.class) {
                setWindowFunction(WindowFunction.none);
            }
        }
        if (properties.getWindowingFunction() != null) {
            setWindowFunction(properties.getWindowingFunction());
        }
        if (properties.getUrl() != null) {
            setFeatureInfoURL(properties.getUrl());
        }

        Map<String, String> attributes = properties.getAttributes();
        if (attributes != null) {
            for (Map.Entry<String, String> entry : attributes.entrySet()) {
                this.setAttributeValue(entry.getKey(), entry.getValue());
            }
        }

        Boolean as = properties.getAutoScale();
        if (as != null) {
            this.autoScale = as;
        }

        if (properties.getFeatureVisibilityWindow() >= 0) {
            setVisibilityWindow(properties.getFeatureVisibilityWindow());
        }
    }

    /**
     * @return the top
     */
    public int getY() {
        return top;
    }

    public void setColorScale(ContinuousColorScale colorScale) {
        this.colorScale = colorScale;
    }

    /**
     * @param top the top to set
     */
    public void setY(int top) {
        this.top = top;
    }

    /**
     * Return the color scale for this track.  Used for heatmaps.
     *
     * @return
     */
    public ContinuousColorScale getColorScale() {

        if (colorScale == null) {

            // Check for a default color scale for this track type in the session.  This is a long deprecated
            // session feature,  but we need to support it for backward compatibility.
            ContinuousColorScale defaultScale = IGV.getInstance().getSession().getColorScale(trackType);
            if (defaultScale == null) {
                // Now check user preferences
                defaultScale = PreferencesManager.getPreferences().getColorScale(trackType);
            }
            if (defaultScale != null) {
                return defaultScale;
            }

            // There is no default,  create one from the data range and track colors
            double min = dataRange == null ? 0 : dataRange.getMinimum();
            double max = dataRange == null ? 10 : dataRange.getMaximum();
            if (min < 0) {
                Color minColor = altColor == null ? oppositeColor(Color.white) : altColor;
                colorScale = new ContinuousColorScale(min, 0, max, minColor, Color.white, getColor());
            } else {
                colorScale = new ContinuousColorScale(min, max, Color.white, getColor());
            }
            colorScale.setNoDataColor(PreferencesManager.getPreferences().getAsColor(Constants.NO_DATA_COLOR));
        }

        return colorScale;
    }

    private Color oppositeColor(Color c) {
        float[] rgb = new float[4];
        c.getRGBComponents(rgb);
        rgb[0] = Math.abs(rgb[0] - 255);
        rgb[1] = Math.abs(rgb[1] - 255);
        rgb[2] = Math.abs(rgb[2] - 255);
        return Color.getHSBColor(rgb[0], rgb[1], rgb[2]);
    }


    public boolean isItemRGB() {
        return itemRGB;
    }

    public boolean isUseScore() {
        return useScore;
    }

    public int getFontSize() {
        return fontSize;
    }

    public void setFontSize(int fontSize) {
        this.fontSize = fontSize;
    }

    public boolean isShowDataRange() {
        return showDataRange;
    }

    public void setShowDataRange(boolean showDataRange) {
        this.showDataRange = showDataRange;
    }


    /**
     * Overriden by subclasses
     *
     * @param e
     * @return
     */
    public Feature getFeatureAtMousePosition(TrackClickEvent e) {
        return null;
    }


    public float logScaleData(float dataY) {

        if (Float.isNaN(dataY)) {
            return dataY;
        }

        // Special case for copy # -- centers data around 2 copies (1 for allele
        // specific) and log normalizes
        if (((getTrackType() == TrackType.COPY_NUMBER) ||
                (getTrackType() == TrackType.ALLELE_SPECIFIC_COPY_NUMBER) ||
                (getTrackType() == TrackType.CNV)) &&
                !isLogNormalized()) {
            double centerValue = (getTrackType() == TrackType.ALLELE_SPECIFIC_COPY_NUMBER)
                    ? 1.0 : 2.0;

            return (float) (Globals.log2(Math.max(Float.MIN_VALUE, dataY) / centerValue));
        } else {
            return dataY;
        }
    }

    public boolean isRegionScoreType(RegionScoreType type) {
        return (getTrackType() == TrackType.GENE_EXPRESSION && type == RegionScoreType.EXPRESSION) ||
                ((getTrackType() == TrackType.COPY_NUMBER || getTrackType() == TrackType.CNV ||
                        getTrackType() == TrackType.ALLELE_SPECIFIC_COPY_NUMBER) &&
                        (type == RegionScoreType.AMPLIFICATION ||
                                type == RegionScoreType.DELETION ||
                                type == RegionScoreType.FLUX)) ||
                (type == RegionScoreType.MUTATION_COUNT) ||
                (type == RegionScoreType.SCORE);
    }

    public void setVisibilityWindow(int i) {
        this.visibilityWindow = i;
    }

    public int getVisibilityWindow() {
        return visibilityWindow;
    }

    /**
     * Override to return a specialized popup menu
     *
     * @return
     */
    public IGVPopupMenu getPopupMenu(final TrackClickEvent te) {
        return null;
    }

    public DisplayMode getDisplayMode() {
        return displayMode;
    }

    public void setDisplayMode(DisplayMode mode) {
        this.displayMode = mode;
    }


    public String getTooltipText(int y) {

        StringBuffer buffer = new StringBuffer();
        buffer.append("<html>" + getName());

        if (resourceLocator != null) {
            Map<String, String> metadata = resourceLocator.getMetadata();
            if (metadata != null && metadata.size() > 0) {
                for (Map.Entry<String, String> entry : metadata.entrySet()) {
                    String value = entry.getValue();
                    if (value != null && value.length() > 0) {
                        buffer.append("<br>" + entry.getKey() + ": " + entry.getValue());
                    }
                }
            }
        }

        if (resourceLocator != null && resourceLocator.getPath() != null) {
            buffer.append("<br>" + this.resourceLocator.getPath());
        }

        return buffer.toString();
    }

    /**
     * Return a value string for the tooltip window at the given location, or null to signal there is no value
     * at that location
     *
     * @param chr
     * @param position
     * @param mouseX
     * @param frame    @return
     */
    public String getValueStringAt(String chr, double position, int mouseX, int mouseY, ReferenceFrame frame) {
        return null;
    }


    public void setWindowFunction(WindowFunction type) {
        // Required method for track interface, ignore
    }

    public WindowFunction getWindowFunction() {
        return null;
    }

    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, String frameName) {
        // Required method for track interface, ignore
        return getRegionScore(chr, start, end, zoom, type, frameName, null);
    }


    /**
     * @param chr
     * @param start
     * @param end
     * @param zoom
     * @param type
     * @param frameName
     * @param tracks
     * @return
     */
    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, String frameName, List<Track> tracks) {
        // Required method for track interface, ignore
        return 0;
    }


    public boolean isLogNormalized() {
        // Required method for track interface, ignore
        return true;
    }

    public boolean isSortable() {
        return sortable;
    }

    public void setSortable(boolean sortable) {
        this.sortable = sortable;
    }

    public boolean isDrawYLine() {
        return drawYLine;
    }

    public float getYLine() {
        return yLine;
    }

    @Override
    public void unload() {
        if (this instanceof IGVEventObserver) {
            IGVEventBus.getInstance().unsubscribe((IGVEventObserver) this);
        }
    }


    public void setRenderer(Renderer renderer) {
        //Here as setter for corresponding getter, subclasses should override
    }

    @Override
    public Renderer getRenderer() {
        return null;
    }

    // Start of Roche-Tessella modification
    public boolean getAutoScale() {
        return this.autoScale;
    }
    // End of Roche-Tessella modification


    public void setShowFeatureNames(boolean b) {
        this.showFeatureNames = b;
    }

    @Override
    public boolean isShowFeatureNames() {
        return showFeatureNames;
    }

    @Override
    public void setTrackLine(String trackLine) {
        this.trackLine = trackLine;
    }

    /**
     * Return "track" line information for exporting features to a file.  Default is null, subclasses may override.
     *
     * @return
     */
    public String getExportTrackLine() {
        return trackLine;
    }

    /**
     * Restore track from XML serialization -- work in progress
     * //        <renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count">
     *
     * @param element
     */

    @Override
    public void marshalXML(Document document, Element element) {

        element.setAttribute("name", name);
        element.setAttribute("attributeKey", attributeKey);
        element.setAttribute("id", id);
        element.setAttribute("fontSize", String.valueOf(fontSize));
        element.setAttribute("visible", String.valueOf(visible));

        if (showFeatureNames != DEFAULT_SHOW_FEATURE_NAMES) {
            element.setAttribute("showFeatureNames", Boolean.toString(showFeatureNames));
        }
        if (color != null) {
            element.setAttribute(SessionAttribute.COLOR, ColorUtilities.colorToString(color));
        }
        if (altColor != null) {
            element.setAttribute(SessionAttribute.ALT_COLOR, ColorUtilities.colorToString(altColor));
        }
        if (visibilityWindow != VISIBILITY_WINDOW) {
            element.setAttribute("featureVisibilityWindow", String.valueOf(visibilityWindow));
        }
        if (displayMode != DEFAULT_DISPLAY_MODE) {
            element.setAttribute(SessionAttribute.DISPLAY_MODE, displayMode.toString());
        }
        if (height != DEFAULT_HEIGHT) {
            element.setAttribute(SessionAttribute.HEIGHT, height.toString());
        }
        if (colorScale != null) {
            //colorScale="ContinuousColorScale;-0.1;-1.5;0.1;1.5;0,153,204;255,255,255;255,0,0"
            element.setAttribute("colorScale", colorScale.asString());
        }
        if (height != DEFAULT_HEIGHT) {
            element.setAttribute("height", String.valueOf(this.height));
        }
        if (showDataRange == false) {
            element.setAttribute("showDataRange", Boolean.toString(showDataRange));
        }
        if (isNumeric()) {
            if (autoscaleGroup != null) {
                element.setAttribute("autoscaleGroup", this.autoscaleGroup);
            }

            element.setAttribute("autoScale", String.valueOf(this.autoScale));

            if (this.getWindowFunction() != null) {
                element.setAttribute("windowFunction", String.valueOf(this.getWindowFunction()));
            }
        }

    }


    @Override
    public void unmarshalXML(Element element, Integer version) {

        if (element.hasAttribute("name")) {
            this.name = element.getAttribute("name");
        }

        if (element.hasAttribute("id")) {
            this.id = element.getAttribute("id");
        }

        if (element.hasAttribute("attributeKey")) {
            this.attributeKey = element.getAttribute("attributeKey");
        } else {
            this.attributeKey = this.name;
        }

        if (element.hasAttribute("displayMode")) {
            try {
                this.displayMode = DisplayMode.valueOf(element.getAttribute("displayMode"));
            } catch (IllegalArgumentException e) {
                log.error("Unrecognized displayMode: " + element.getAttribute("displayMode"));
                this.displayMode = DisplayMode.COLLAPSED;
            }
        }

        if (element.hasAttribute("color")) {
            try {
                Color c = ColorUtilities.stringToColor(element.getAttribute("color"));
                this.color = c;
            } catch (Exception e) {
                log.error("Unrecognized color: " + element.getAttribute("color"));
            }
        }

        if (element.hasAttribute("altColor")) {
            try {
                Color c = ColorUtilities.stringToColor(element.getAttribute("altColor"));
                this.altColor = c;
            } catch (Exception e) {
                log.error("Unrecognized altColor: " + element.getAttribute("altColor"));
            }
        }

        if (element.hasAttribute("colorScale")) {
            try {
                this.colorScale = (ContinuousColorScale) ColorScaleFactory.getScaleFromString(element.getAttribute("colorScale"));
            } catch (Exception e) {
                log.error("Unrecognized colorScale: " + element.getAttribute("colorScale"));
            }
        }

        if (element.hasAttribute("visible")) {
            try {
                this.setVisible(Boolean.parseBoolean(element.getAttribute("visible")));
            } catch (Exception e) {
                log.error("Unrecognized visisbilty: " + element.getAttribute("visible"));
            }
        }

        if (element.hasAttribute("autoScale")) {
            try {
                this.autoScale = Boolean.valueOf(element.getAttribute("autoScale"));
            } catch (Exception e) {
                log.error("Unrecognized autoScale: " + element.getAttribute("autoScale"));
            }
        }

        if (element.hasAttribute("autoscaleGroup")) {
            String autoscaleGroup = element.getAttribute("autoscaleGroup");
            this.setAttributeValue(AttributeManager.GROUP_AUTOSCALE, "" + autoscaleGroup);
        }

        if (element.hasAttribute("showDataRange")) {
            try {
                this.showDataRange = Boolean.valueOf(element.getAttribute("showDataRange"));
            } catch (Exception e) {
                log.error("Unrecognized showDataRange: " + element.getAttribute("showDataRange"));
            }
        }

        if (element.hasAttribute("featureVisibilityWindow")) {
            try {
                this.visibilityWindow = Integer.parseInt(element.getAttribute("featureVisibilityWindow"));
            } catch (NumberFormatException e) {
                log.error("Unrecognized featureVisibilityWindow: " + element.getAttribute("featureVisibilityWindow"));
            }
        }

        if (element.hasAttribute("showFeatureNames")) {
            try {
                this.showFeatureNames = Boolean.valueOf(element.getAttribute("showFeatureNames"));
            } catch (Exception e) {
                log.error("Unrecognized showDataRange: " + element.getAttribute("showFeatureNames"));
            }

        }

        if (element.hasAttribute("fontSize")) {
            try {
                this.fontSize = Integer.parseInt(element.getAttribute("fontSize"));
            } catch (NumberFormatException e) {
                log.error("Unrecognized fontSize: " + element.getAttribute("fontSize"));
            }
        }

        if (element.hasAttribute("height")) {
            try {
                this.height = Integer.parseInt(element.getAttribute("height"));
            } catch (NumberFormatException e) {
                log.error("Unrecognized height: " + element.getAttribute("height"));
            }
        }

        if (element.hasAttribute("windowFunction")) {
            try {
                this.setWindowFunction(WindowFunction.valueOf(element.getAttribute("windowFunction")));
            } catch (IllegalArgumentException e) {
                log.error("Unknown windowFunction: " + element.getAttribute("windowFunction"), e);
            }
        }

        // Set DataRange -- legacy (pre V3 sessions)
        if (version <= 3 && element.hasAttribute(SessionAttribute.SCALE)) {
            String scale = element.getAttribute(SessionAttribute.SCALE);
            try {
                String[] axis = scale.split(",");
                float minimum = Float.parseFloat(axis[0]);
                float baseline = Float.parseFloat(axis[1]);
                float maximum = Float.parseFloat(axis[2]);
                setDataRange(new DataRange(minimum, baseline, maximum));
            } catch (NumberFormatException e) {
                log.error("Unrecognized dataRange: " + element.getAttribute("scale"));
            }
        }

        NodeList nodeList = element.getElementsByTagName("DataRange");
        if (nodeList != null && nodeList.getLength() > 0) {
            Element dataRangeElement = (Element) nodeList.item(0);
            try {
                this.dataRange = new DataRange(dataRangeElement, version);
            } catch (Exception e) {
                log.error("Unrecognized DataRange");
            }
        }
    }

    public String getSampleId() {
        return sampleId;
    }

    public Color getFeatureColor(IGVFeature feature) {

        // Set color used to draw the feature
        Color featureColor = null;

        // If an alt color is explicitly set use it for negative strand features;
        if (feature.getStrand() == Strand.NEGATIVE) {
            featureColor = getAltColor();
        }

        // If color is explicitly set use it
        if (featureColor == null) {
            featureColor = getColor();
        }

        // No explicitly set color, try the feature itself
        if (featureColor == null) {
            featureColor = feature.getColor();
        }

        // If still no color use defaults
        if (featureColor == null) {
            if (getTrackType() == TrackType.CNV) {
                featureColor = feature.getName().equals("gain") ? Globals.DULL_RED : Globals.DULL_BLUE;
            } else {
                featureColor = getDefaultColor();
            }
        }

        if (useScore) {
            float score = feature.getScore();
            float alpha = 1;
            if (!Float.isNaN(score)) {
                float binWidth = (viewLimitMax - viewLimitMin) / 9;
                int binNumber = (int) ((score - viewLimitMin) / binWidth);
                alpha = Math.min(1.0f, 0.2f + (binNumber * 0.8f) / 9);
            }
            featureColor = ColorUtilities.getCompositeColor(featureColor, alpha);
        }

        return featureColor;
    }

}
