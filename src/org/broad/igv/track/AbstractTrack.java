/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.track;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.gwas.GWASTrack;
import org.broad.igv.renderer.*;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.sam.CoverageTrack;
import org.broad.igv.session.IGVSessionReader;
import org.broad.igv.session.SessionXmlAdapters;
import org.broad.igv.session.SubtlyImportant;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.TooltipTextFrame;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.panel.AttributeHeaderPanel;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.MouseableRegion;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.Utilities;
import htsjdk.tribble.Feature;
import org.w3c.dom.Node;

import javax.xml.bind.JAXBException;
import javax.xml.bind.annotation.*;
import javax.xml.bind.annotation.adapters.XmlJavaTypeAdapter;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 */
@XmlType(factoryClass = IGVSessionReader.class, factoryMethod = "getNextTrack")
@XmlAccessorType(XmlAccessType.NONE)
@XmlSeeAlso({CoverageTrack.class, AlignmentTrack.class, DataSourceTrack.class, GWASTrack.class, FeatureTrack.class, MergedTracks.class})
public abstract class AbstractTrack implements Track {

    private static Logger log = Logger.getLogger(AbstractTrack.class);

    /**
     * Classes which we have tried to marshal/unmarshal
     * and have failed. Since we just use exception catching (slow),
     * we don't want to repeat failures
     */
    public static final Set<Class> knownUnknownTrackClasses = new HashSet<Class>();
    public static final Class defaultTrackClass = AbstractTrack.class;

    /**
     * Set default renderer classes by track type.
     */
    private static Class defaultRendererClass = BarChartRenderer.class;
    private static Map<TrackType, Class> defaultRendererMap = new HashMap();

    static {
        defaultRendererMap.put(TrackType.RNAI, HeatmapRenderer.class);
        defaultRendererMap.put(TrackType.COPY_NUMBER, HeatmapRenderer.class);
        defaultRendererMap.put(TrackType.CNV, HeatmapRenderer.class);
        defaultRendererMap.put(TrackType.ALLELE_SPECIFIC_COPY_NUMBER, HeatmapRenderer.class);
        defaultRendererMap.put(TrackType.GENE_EXPRESSION, HeatmapRenderer.class);
        defaultRendererMap.put(TrackType.DNA_METHYLATION, HeatmapRenderer.class);
        defaultRendererMap.put(TrackType.LOH, HeatmapRenderer.class);
        defaultRendererMap.put(TrackType.OTHER, BarChartRenderer.class);
        defaultRendererMap.put(TrackType.CHIP_CHIP, HeatmapRenderer.class);
    }


    @XmlAttribute
    private String id;
    @XmlAttribute
    private String name;
    private String url;
    private boolean itemRGB = true;

    private boolean useScore;
    private float viewLimitMin = Float.NaN;     // From UCSC track line
    private float viewLimitMax = Float.NaN;  // From UCSC track line


    @XmlAttribute
    protected int fontSize = PreferenceManager.getInstance().getAsInt(PreferenceManager.DEFAULT_FONT_SIZE);
    private boolean showDataRange = true;
    private String sampleId;

    private ResourceLocator resourceLocator;

    private int top;
    protected int minimumHeight = -1;
    protected int maximumHeight = 1000;

    private TrackType trackType = TrackType.OTHER;

    private boolean selected = false;

    @XmlAttribute
    private boolean visible = true;
    @XmlAttribute
    private boolean sortable = true;

    boolean overlaid;

    boolean drawYLine = false;
    float yLine = 0;

    // Map to store attributes specific to this track.  Attributes shared by multiple
    private Map<String, String> attributes = new HashMap();

    // Scale for heatmaps
    @XmlJavaTypeAdapter(SessionXmlAdapters.ContinuousColorScale.class)
    @XmlAttribute
    private ContinuousColorScale colorScale;

    //TODO Only write it out for applicable tracks
    //Not applicable to all tracks.
    @XmlAttribute
    protected boolean autoScale;

    @XmlJavaTypeAdapter(SessionXmlAdapters.Color.class)
    @XmlAttribute(name = "color")
    private Color posColor = Color.blue.darker();

    @XmlJavaTypeAdapter(SessionXmlAdapters.Color.class)
    @XmlAttribute
    private Color altColor = posColor;

    @XmlAttribute(name = "featureVisibilityWindow")
    protected int visibilityWindow = -1;
    @XmlAttribute
    private DisplayMode displayMode = DisplayMode.COLLAPSED;

    @XmlJavaTypeAdapter(SessionXmlAdapters.Height.class)
    @XmlAttribute
    protected Integer height = -1;

    @XmlElement(name = "DataRange")
    protected DataRange dataRange;

    @SubtlyImportant
    private AbstractTrack() {
    }

    public AbstractTrack(
            ResourceLocator dataResourceLocator,
            String id,
            String name) {
        this.resourceLocator = dataResourceLocator;
        this.id = id;
        this.name = name;
        init();
    }

    public AbstractTrack(ResourceLocator dataResourceLocator, String id) {
        this(dataResourceLocator, id, dataResourceLocator.getTrackName());
    }

    public AbstractTrack(ResourceLocator locator) {
        this(locator, locator != null ? locator.getPath() : null, locator != null ? locator.getTrackName() : null);
    }

    public AbstractTrack(String id) {
        this(null, id, id);
    }


    public AbstractTrack(String id, String name) {
        this(null, id, name);
    }

    private void init() {
        showDataRange = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.CHART_SHOW_DATA_RANGE);
        if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.EXPAND_FEAUTRE_TRACKS)) {
            displayMode = DisplayMode.EXPANDED;
        }
    }

    public void setRendererClass(Class rc) {
        // Ignore by default
    }

    public String getUrl() {
        return url;
    }

    public void setUrl(String url) {
        this.url = url;
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

    private String getDisplayName() {

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

    @Override
    public void load(RenderContext context) {
        // No-op, to be overriden by subclasses
    }

    @Override
    public boolean isFilterable() {
        return true;   // True by default
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


    protected int lastRenderY = -1;
    /**
     * Render the components which should only be rendered once
     * for a give y-value. We sometimes render tracks piecewise, but may
     * only want to show the scale once at the left, for instance
     * @param context
     * @param rect
     */

    /**
     * Return whether we have rendered anything at this Y-coordinate already
     * @param rect
     * @return
     */
    protected boolean isRepeatY(Rectangle rect){
        return rect.y == lastRenderY;
    }

    public void resetLastY(){
        this.lastRenderY = -1;
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


    public Color getColor() {
        return posColor;
    }

    public Color getAltColor() {
        return altColor;

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
        attributes.put(key, value);
        AttributeManager.getInstance().addAttribute(getSample(), name, value);
    }


    /**
     * Return the attribute value.  Attribute lookup occurs in the following order, if all fail null is returned.
     * <p/>
     * (1) the track attribute table
     * (2) by sampleId, as set in the Resource element of a session or load-from-server menu
     * (3) by track name, the visible display name
     * (4) by full path to the file associated with this track
     *
     * @param attributeName
     * @return
     */
    public String getAttributeValue(String attributeName) {
        String key = attributeName.toUpperCase();
        String value = attributes.get(key);
        final AttributeManager attributeManager = AttributeManager.getInstance();
        if (value == null && getSample() != null) {
            value = attributeManager.getAttribute(getSample(), key);
        }
        if (value == null) {
            value = attributeManager.getAttribute(getName(), key);
        }
        if (value == null && getResourceLocator() != null && getResourceLocator().getPath() != null) {
            value = attributeManager.getAttribute(getResourceLocator().getPath(), key);
        }
        return value;
    }

    public String getSample() {

        if (sampleId != null) {
            return sampleId;
        }
//        String sample = AttributeManager.getInstance().getSampleFor(getName());
//        return sample != null ? sample : getName();

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
        if (XYPlotRenderer.class.isAssignableFrom(getDefaultRendererClass())) {
            return PreferenceManager.getInstance().getAsInt(PreferenceManager.CHART_TRACK_HEIGHT_KEY);
        } else {
            return PreferenceManager.getInstance().getAsInt(PreferenceManager.TRACK_HEIGHT_KEY);
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

    public void setMaximumHeight(int maximumHeight) {
        this.maximumHeight = maximumHeight;
    }


    /**
     * Return the actual minimum height if one has been set, otherwise get the default for the current renderer.
     *
     * @return
     */
    public int getMinimumHeight() {
        return minimumHeight < 0 ? getDefaultMinimumHeight() : minimumHeight;
    }

    public int getMaximumHeight() {
        return maximumHeight;
    }

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
                    return PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_ORPHANED_MUTATIONS);
                }
            }
        }
        return visible;
    }

    public void setColor(Color color) {
        this.posColor = color;
    }


    public void setAltColor(Color color) {
        altColor = color;
    }


    public void setVisible(boolean isVisible) {
        this.visible = isVisible;
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
            this.height = Math.min(Math.max(getMinimumHeight(), preferredHeight), getMaximumHeight());
        }
    }

    public int getHeight() {
        return (height < 0) ? getDefaultHeight() : height;
    }

    public boolean hasDataRange() {
        return dataRange != null;
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


    protected Class getDefaultRendererClass() {
        Class def = defaultRendererMap.get(getTrackType());
        return (def == null) ? defaultRendererClass : def;
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
        String popupText = getValueStringAt(frame.getChrName(), e.getChromosomePosition(), e.getMouseEvent().getY(), frame);

        if (popupText != null) {

            final TooltipTextFrame tf = new TooltipTextFrame(getName(), popupText);
            Point p = me.getComponent().getLocationOnScreen();
            tf.setLocation(Math.max(0, p.x + me.getX() - 150), Math.max(0, p.y + me.getY() - 150));

            UIUtilities.invokeOnEventThread(new Runnable() {
                public void run() {
                    tf.setVisible(true);
                }
            });
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
        this.itemRGB = properties.isItemRGB();
        this.useScore = properties.isUseScore();
        this.viewLimitMin = properties.getMinValue();
        this.viewLimitMax = properties.getMaxValue();
        this.yLine = properties.getyLine();
        this.drawYLine = properties.isDrawYLine();
        this.sortable = properties.isSortable();


        // If view limits are explicitly set turn off autoscale
        if (!Float.isNaN(viewLimitMin) && !Float.isNaN(viewLimitMax)) {
            this.setAutoScale(false);
        } else {
            this.setAutoScale(properties.isAutoScale());
        }

        // Color scale properties
        if (!properties.isAutoScale()) {

            float min = properties.getMinValue();
            float max = properties.getMaxValue();

            float mid = properties.getMidValue();
            if (Float.isNaN(mid)) {
                if (min >= 0) {
                    mid = Math.max(min, 0);
                } else {
                    mid = Math.min(max, 0);
                }
            }


            DataRange dr = new DataRange(min, mid, max);
            setDataRange(dr);

            if (properties.isLogScale()) {
                dr.setType(DataRange.Type.LOG);
            }

            // If the user has explicity set a data range and colors apply to heatmap as well
            Color maxColor = properties.getColor();
            Color minColor = properties.getAltColor();
            if (maxColor != null && minColor != null) {

                float tmp = properties.getNeutralFromValue();
                float neutralFrom = Float.isNaN(tmp) ? mid : tmp;
                tmp = properties.getNeutralToValue();
                float neutralTo = Float.isNaN(tmp) ? mid : tmp;

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
        if (properties.getHeight() > 0) {
            setHeight(properties.getHeight());
        }
        if (properties.getMinHeight() > 0) {
            setMinimumHeight(properties.getMinHeight());
        }
        if (properties.getRendererClass() != null) {
            setRendererClass(properties.getRendererClass());
        }
        if (properties.getWindowingFunction() != null) {
            setWindowFunction(properties.getWindowingFunction());
        }
        if (properties.getUrl() != null) {
            setUrl(properties.getUrl());
        }

        Map<String, String> attributes = properties.getAttributes();
        if (attributes != null) {
            for (Map.Entry<String, String> entry : attributes.entrySet()) {
                this.setAttributeValue(entry.getKey(), entry.getValue());
            }
        }

        // Start of Roche-Tessella modification
        this.autoScale = properties.getAutoScale();
        // End of Roche-Tessella modification

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
     * Return the color scale for this track.  If a specific scale exists for this data type
     * use that.  Otherwise create one using the track color and data range.
     *
     * @return
     */
    public ContinuousColorScale getColorScale() {
        if (colorScale == null) {

            if (IGV.hasInstance()) {
                ContinuousColorScale defaultScale = IGV.getInstance().getSession().getColorScale(trackType);
                if (defaultScale != null) {
                    return defaultScale;
                }
            }

            double min = dataRange == null ? 0 : dataRange.getMinimum();
            double max = dataRange == null ? 10 : dataRange.getMaximum();
            Color c = getColor();
            Color minColor = Color.white;
            if (min < 0) {
                minColor = altColor == null ? oppositeColor(minColor) : altColor;
                colorScale = new ContinuousColorScale(min, 0, max, minColor, Color.white, c);
            } else {
                colorScale = new ContinuousColorScale(min, max, minColor, c);
            }
            colorScale.setNoDataColor(UIConstants.NO_DATA_COLOR);
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

    /**
     * Return the current state of this object as map of key-value pairs.  Used to store session state.
     * Only those attributes not already annotated in AbstractTrack need to be included here
     *
     * @return
     * @see #restorePersistentState
     */
    public Map<String, String> getPersistentState() {
        return new HashMap<String, String>();
    }


    /**
     * Restore from XML node. Default implementation just turns attributes
     * into a map
     *
     * @param node
     */
    public void restorePersistentState(Node node, int version) throws JAXBException {
        Map<String, String> attributes = Utilities.getAttributes(node);
        restorePersistentState(attributes, version);
    }

    /**
     * Restore attributes from track tag, no children
     * Only those attributes not unmarshalled (meaning not part of AbstractTrack)
     * need to be restored
     *
     * @param attributes
     * @see #getPersistentState
     */
    public void restorePersistentState(Map<String, String> attributes, int version) {

        String displayName = attributes.get(IGVSessionReader.SessionAttribute.DISPLAY_NAME.getText());
        String name = attributes.get(IGVSessionReader.SessionAttribute.NAME.getText());

        if (name != null && name.length() > 0) {
            setName(name);
        } else if (displayName != null && displayName.length() > 0) {
            setName(displayName);
        }

        // Set DataRange -- legacy (pre V3 sessions)
        if(version <= 3){
            String scale = attributes.get(IGVSessionReader.SessionAttribute.SCALE.getText());
            if (scale != null) {
                String[] axis = scale.split(",");
                float minimum = Float.parseFloat(axis[0]);
                float baseline = Float.parseFloat(axis[1]);
                float maximum = Float.parseFloat(axis[2]);
                setDataRange(new DataRange(minimum, baseline, maximum));
            }
        }

    }

    public boolean isItemRGB() {
        return itemRGB;
    }

    public boolean isUseScore() {
        return useScore;
    }

    public float getViewLimitMin() {
        return viewLimitMin;
    }

    public float getViewLimitMax() {
        return viewLimitMax;
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

    /**
     * Special normalization function for linear (non logged) copy number data
     *
     * @param value
     * @param norm
     * @return
     */
    public static float getLogNormalizedValue(float value, double norm) {
        if (norm == 0) {
            return Float.NaN;
        } else {
            return (float) (Math.log(Math.max(Float.MIN_VALUE, value) / norm) / Globals.log2);
        }
    }

    public float logScaleData(float dataY) {

        // Special case for copy # -- centers data around 2 copies (1 for allele
        // specific) and log normalizes
        if (((getTrackType() == TrackType.COPY_NUMBER) ||
                (getTrackType() == TrackType.ALLELE_SPECIFIC_COPY_NUMBER) ||
                (getTrackType() == TrackType.CNV)) &&
                !isLogNormalized()) {
            double centerValue = (getTrackType() == TrackType.ALLELE_SPECIFIC_COPY_NUMBER)
                    ? 1.0 : 2.0;

            dataY = getLogNormalizedValue(dataY, centerValue);
        }


        return dataY;
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


    public String getNameValueString(int y) {
        return getName();
    }

    /**
     * Return a value string for the tooltip window at the given location, or null to signal there is no value
     * at that location
     *
     * @param chr
     * @param position
     * @param y
     * @param frame
     * @return
     */
    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {
        return null;
    }


    public void setWindowFunction(WindowFunction type) {
        // Required method for track interface, ignore
    }

    @XmlAttribute(name = "windowFunction")
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
    public void dispose() {
        // Default is to do nothing.  Override in subclasses

    }


    protected void setRenderer(Renderer renderer) {
        //Here as setter for corresponding getter, subclasses should override
    }

    @XmlJavaTypeAdapter(SessionXmlAdapters.Renderer.class)
    @XmlAttribute(name = "renderer")
    @Override
    public Renderer getRenderer() {
        return null;
    }

    // Start of Roche-Tessella modification
    public boolean getAutoScale()    {
        return this.autoScale;
    }
    // End of Roche-Tessella modification
}
