/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */
package org.broad.igv.track;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.session.RendererFactory;
import org.broad.igv.ui.IGVMainFrame;
import org.broad.tribble.Feature;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.renderer.*;
import org.broad.igv.session.SessionReader;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.util.ColorUtilities;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 */
public abstract class AbstractTrack implements Track {

    private static Logger log = Logger.getLogger(AbstractTrack.class);

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


    private String id;
    private String id_142;    // 1.4.2 name, used to restore old sessions.
    private String name;
    private Color midColor;
    private String url;
    private boolean itemRGB = true;
    private boolean useScore;
    private float viewLimitMin;
    private float viewLimitMax;
    private int fontSize = 9;
    private boolean showDataRange = true;
    private String sampleId;

    /**
     * The source file from which this track was created.  Defaults to an
     * empty string, which indicates that the track was not derived from a
     * data file.
     */
    private String sourceFile = "";
    /**
     * Top (Y) position
     */
    private int top;
    /**
     * Minimum allowable height for this track
     */
    private int minimumHeight = 1;
    /**
     * The file from which this track originated.
     */
    private ResourceLocator resourceLocator;

    private TrackType trackType = TrackType.OTHER;

    /**
     * Continuous posColor scale for this track
     *
     * @return
     */
    private boolean isSelected = false;
    private boolean visible = true;
    boolean overlayVisible;
    private int preferredHeight;
    private boolean isDraggable = true;

    /**
     * Map to store attributes specific to this track.  Attributes shared by multiple
     * tracks are stored in AttributeManager.
     */
    private Map<String, String> attributes = new HashMap();
    /**
     * Scale for heatmaps
     */
    private ContinuousColorScale colorScale;
    /**
     * Chart parameters
     */
    private Color posColor;
    private Color altColor;
    private DataRange dataRange;
    protected int visibilityWindow = -1;


    public AbstractTrack(
            ResourceLocator dataResourceLocator,
            String id,
            String name) {
        this.resourceLocator = dataResourceLocator;
        this.id = id;
        this.name = name;
        init();
    }

    public AbstractTrack(ResourceLocator dataResourceLocator) {
        this.resourceLocator = dataResourceLocator;
        this.id = dataResourceLocator.getPath();

        String drName = dataResourceLocator.getName();
        this.name = drName != null ? drName : dataResourceLocator.getFileName();


        init();
    }


    public AbstractTrack(String id) {
        this.name = id;
        this.id = id;
        init();

    }


    public AbstractTrack(String id, String name) {
        this.name = name;
        this.id = id;
        init();
    }

    private void init() {
        overlayVisible = IGVMainFrame.getInstance().getSession().getDisplayOverlayTracks();
        showDataRange = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.CHART_SHOW_DATA_RANGE);
        this.id_142 = String.valueOf(name == null ? resourceLocator.getPath() : name);

    }


    public String getUrl() {
        return url;
    }

    public void setUrl(String url) {
        this.url = url;
    }


    public String getId() {
        return id;
    }

    /**
     * Method description
     *
     * @param name
     */
    public void setName(String name) {
        this.name = name;
    }

    // TODO Provided for session saving.  Rename/refactor later


    /**
     * Method description
     *
     * @return
     */
    public String getName() {

        return name;
    }

    private String getDisplayName() {

        String sampleKey = IGVMainFrame.getInstance().getSession().getTrackAttributeName();
        if (sampleKey != null && sampleKey.trim().length() > 0) {
            String name = getAttributeValue(sampleKey.trim());
            if (name != null) {
                return name;
            }
        }
        return getName();
    }


    public String getSampleId() {
        if (sampleId != null) {
            return sampleId;
        } else {
            return getName();
        }
    }

    public void setSampleId(String sampleId) {
        this.sampleId = sampleId;
    }

    /**
     * Render a border.
     *
     * @param context
     * @param rect
     */
    public void renderBorder(RenderContext context, Rectangle rect) {
        Renderer renderer = this.getRenderer();
        if (renderer != null) {
            this.getRenderer().renderBorder(this, context, rect);
        }
    }

    /**
     * Render a Y axis
     *
     * @param context
     * @param rect
     */
    public void renderAxis(RenderContext context, Rectangle rect) {
        Renderer renderer = this.getRenderer();
        if (renderer != null) {
            this.getRenderer().renderAxis(this, context, rect);
        }
    }

    public void renderName(Graphics2D graphics, Rectangle trackRectangle, Rectangle visibleRectangle) {

        Rectangle rect = getDisplayableRect(trackRectangle, visibleRectangle);

        String trackName = getDisplayName();
        if ((trackName != null)) {
            if (isSelected()) {
                graphics.setBackground(Color.LIGHT_GRAY);
            } else {
                graphics.setBackground(Color.WHITE);
            }
            Graphics2D g2D = graphics; //(Graphics2D) graphics.create();
            g2D.clearRect(rect.x, rect.y, rect.width, rect.height);
            if (rect.getHeight() > 3) {


                // Calculate fontsize
                int gap = Math.min(4, rect.height / 3);
                int fontSize = Math.min(12, rect.height - gap);

                Font font = FontManager.getScalableFont(fontSize);
                g2D.setFont(font);
                FontManager.applyScalableTextHints(g2D);

                GraphicUtils.drawWrappedText(trackName, rect, g2D, false);

                //g2D.dispose();
            } else {

                /*Font font = FontManager.getScalableFont(2);
                graphics.setFont(font);
                FontMetrics fm = graphics.getFontMetrics();
                int nameWidth = Math.min(rect.width - 2, fm.stringWidth(trackName));

                graphics.setColor(Color.GRAY);
                int y = (rect.y + rect.height / 2);
                int start = rect.x + 1;
                while (start < nameWidth) {
                    graphics.drawLine(start, y, start + 1, y);
                    start += 2 + (int) (2 * Math.random() * 4);
                }
                */

            }
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
        if (posColor != null) {
            return posColor;
        }
        return getRenderer() == null ? null : getRenderer().getDefaultColor();
    }

    public Color getAltColor() {
        return altColor;

    }

    private Color getColorAttribute(String key, Color defaultValue) {
        String rgb = this.getAttributeValue(key);
        if (rgb != null) {
            try {
                return ColorUtilities.convertRGBStringToColor(rgb.replace("\"", ""));

            } catch (Exception exception) {
                log.info("Invalid color string " + rgb + " for track: " + getName());
            }
        }
        return defaultValue;
    }


    public ResourceLocator getResourceLocator() {
        return resourceLocator;
    }

    /**
     * Add an attribute to this track and register the key with the attribute panel.
     * <p/>
     * Note:  Attribute keys are case insensitive.  Currently this is implemented
     * by forcing all keys to upper case
     *
     * @param key
     * @param value
     */
    public void setAttributeValue(String key, String value) {
        String uppercaseKey = key.toUpperCase();
        attributes.put(uppercaseKey, value);
        AttributeManager.getInstance().addAttributeKey(uppercaseKey);
    }


    public String getAttributeValue(String attributeKey) {
        String value = attributes.get(attributeKey);
        if (value == null) {
            value = AttributeManager.getInstance().getAttribute(getName(), attributeKey);
        }
        return value;
    }


    public int getLevelNumber(int y) {
        return (preferredHeight == 0) ? 0 : y / preferredHeight;
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


    public void setMinimumHeight(int minimumHeight) {
        this.minimumHeight = minimumHeight;
    }


    public int getMinimumHeight() {
        return minimumHeight;
    }


    public void setTrackType(TrackType type) {
        this.trackType = type;
    }


    public TrackType getTrackType() {
        return trackType;
    }

    public boolean isVisible() {
        return visible && ((getTrackType() != UIConstants.overlayTrackType) || (this.overlayVisible == true));
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


    public void setOverlayVisible(boolean bool) {
        this.overlayVisible = bool;
    }


    public void updateState() {
    }


    public void setSelected(boolean selected) {
        isSelected = selected;
    }


    public boolean isSelected() {
        return isSelected;
    }


    public void setIsDraggable(boolean value) {
        isDraggable = value;
    }


    public boolean isDraggable() {
        return isDraggable;
    }

    int height = -1;


    public void setHeight(int height) {
        this.height = Math.max(minimumHeight, height);
    }


    public int getHeight() {
        return (height < 0) ? getDefaultHeight() : height;
    }

    public int getPreferredHeight() {
        return getHeight();
    }


    public DataRange getDataRange() {
        if (dataRange == null) {
            // Use the color scale if htere is one
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

    // Assumes features are sorted by start position

    protected Feature getFeatureAt(double position, double minWidth, List<? extends LocusScore> features) {

        return FeatureUtils.getFeatureAt(position, minWidth, features);
    }


    /**
     * Refresh the underlying data for the track.  Default implementation does nothing.
     * Subclasses can override
     *
     * @param timestamp
     */
    public void refreshData(long timestamp) {
    }


    public String getSourceFile() {
        return sourceFile;
    }

    public void setSourceFile(String sourceFile) {
        this.sourceFile = sourceFile;
    }

    boolean expanded = false;

    public void setExpanded(boolean expanded) {
        this.expanded = expanded;
    }

    public boolean isExpanded() {
        return expanded;
    }


    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return new ArrayList();

    }

    public Color getMidColor() {
        return midColor;
    }

    public void setMidColor(Color midColor) {
        this.midColor = midColor;
    }


    public boolean handleClick(TrackClickEvent e) {
        return false;
    }


    public void setTrackProperties(TrackProperties trackProperties) {

        this.itemRGB = trackProperties.isItemRGB();
        this.useScore = trackProperties.isUseScore();
        this.viewLimitMin = trackProperties.getMinValue();
        this.viewLimitMax = trackProperties.getMaxValue();

        // If view limits are explicitly set turn off autoscale
        // TODO -- get rid of this ugly instance of and casting business
        if (!Float.isNaN(viewLimitMin) && !Float.isNaN(viewLimitMax) && (this instanceof DataTrack)) {
            ((DataTrack) this).setAutoscale(false);
        }

        // Color scale properties
        if (!trackProperties.isAutoScale() &&
                !Float.isNaN(trackProperties.getMinValue()) &&
                !Float.isNaN(trackProperties.getMaxValue())) {

            float min = trackProperties.getMinValue();
            float max = trackProperties.getMaxValue();

            float mid = trackProperties.getMidValue();
            if (Float.isNaN(mid)) {
                if (min >= 0) {
                    mid = Math.max(min, 0);
                } else {
                    mid = Math.min(max, 0);
                }
            }


            DataRange dr = new DataRange(min, mid, max);
            setDataRange(dr);

            if (trackProperties.isLogScale()) {
                dr.setType(DataRange.Type.LOG);
            }

            // If the user has explicity set a data range and colors apply to heatmap as well
            Color maxColor = trackProperties.getColor();
            Color minColor = trackProperties.getAltColor();
            if (maxColor != null && minColor != null) {

                float tmp = trackProperties.getNeutralFromValue();
                float neutralFrom = Float.isNaN(tmp) ? mid : tmp;
                tmp = trackProperties.getNeutralToValue();
                float neutralTo = Float.isNaN(tmp) ? mid : tmp;

                Color midColor = trackProperties.getMidColor();
                if (midColor == null) {
                    midColor = Color.white;
                }
                colorScale = new ContinuousColorScale(neutralFrom, min, neutralTo, max, minColor, midColor, maxColor);
            }


        }


        if (trackProperties.getName() != null) {
            name = trackProperties.getName();
        }
        if (trackProperties.getColor() != null) {
            setColor(trackProperties.getColor());
        }
        if (trackProperties.getAltColor() != null) {
            setAltColor(trackProperties.getAltColor());
        }
        if (trackProperties.getMidColor() != null) {
            setMidColor(trackProperties.getMidColor());
        }
        if (trackProperties.getHeight() > 0) {
            setHeight(trackProperties.getHeight());
        }
        if (trackProperties.getMinHeight() > 0) {
            setMinimumHeight(trackProperties.getMinHeight());
        }
        if (trackProperties.getRendererClass() != null) {
            setRendererClass(trackProperties.getRendererClass());
        }
        if (trackProperties.getWindowingFunction() != null) {
            setStatType(trackProperties.getWindowingFunction());
        }
        if (trackProperties.getUrl() != null) {
            setUrl(trackProperties.getUrl());
        }

    }

    /**
     * @return the top
     */
    public int getTop() {
        return top;
    }

    public void setColorScale(ContinuousColorScale colorScale) {
        this.colorScale = colorScale;
    }

    /**
     * @param top the top to set
     */
    public void setTop(int top) {
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

            ContinuousColorScale defaultScale = IGVMainFrame.getInstance().getSession().getColorScale(trackType);
            if (defaultScale != null) {
                return defaultScale;
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


    public Map<String, String> getPersistentState() {

        LinkedHashMap<String, String> attributes = new LinkedHashMap();

        // Color scale
        if (colorScale != null && !colorScale.isDefault()) {
            attributes.put(SessionReader.SessionAttribute.COLOR_SCALE.getText(), colorScale.asString());
        }

        attributes.put("showDataRange", String.valueOf(showDataRange));

        // Visibility
        attributes.put(SessionReader.SessionAttribute.VISIBLE.getText(), String.valueOf(isVisible()));

        // expanded
        attributes.put(SessionReader.SessionAttribute.EXPAND.getText(), String.valueOf(isExpanded()));

        // height
        int height = getHeight();
        if ((this instanceof FeatureTrack) && isExpanded()) {
            height = height / ((FeatureTrack) this).getNumberOfFeatureLevels();
        }
        String value = Integer.toString(height);
        attributes.put(SessionReader.SessionAttribute.HEIGHT.getText(), value);


        if (name != null) {
            attributes.put(SessionReader.SessionAttribute.NAME.getText(), name);

        }

        // color
        Color color = getColor();
        if (color != null) {
            StringBuffer stringBuffer = new StringBuffer();
            stringBuffer.append(color.getRed());
            stringBuffer.append(",");
            stringBuffer.append(color.getGreen());
            stringBuffer.append(",");
            stringBuffer.append(color.getBlue());
            attributes.put(SessionReader.SessionAttribute.COLOR.getText(), stringBuffer.toString());
        }

        // renderer
        Renderer renderer = getRenderer();
        if (renderer != null) {
            RendererFactory.RendererType type = RendererFactory.getRenderType(renderer);
            if (type != null) {
                attributes.put(SessionReader.SessionAttribute.RENDERER.getText(), type.name());
            }
        }

        // window function
        WindowFunction wf = getWindowFunction();
        if (wf != null) {
            attributes.put(SessionReader.SessionAttribute.WINDOW_FUNCTION.getText(), wf.name());
        }


        if (this instanceof FeatureTrack) {
            boolean expand = this.isExpanded();
            attributes.put(SessionReader.SessionAttribute.EXPAND.getText(), String.valueOf(expand));
        }

        attributes.put("fontSize", String.valueOf(getFontSize()));

        return attributes;
    }


    public void restorePersistentState(Map<String, String> attributes) {

        String displayName = attributes.get(SessionReader.SessionAttribute.DISPLAY_NAME.getText());
        String name = attributes.get(SessionReader.SessionAttribute.NAME.getText());

        String isVisible = attributes.get(SessionReader.SessionAttribute.VISIBLE.getText());
        String height = attributes.get(SessionReader.SessionAttribute.HEIGHT.getText());
        String color = attributes.get(SessionReader.SessionAttribute.COLOR.getText());
        String rendererType = attributes.get(SessionReader.SessionAttribute.RENDERER.getText());
        String windowFunction = attributes.get(SessionReader.SessionAttribute.WINDOW_FUNCTION.getText());
        String scale = attributes.get(SessionReader.SessionAttribute.SCALE.getText());
        String isExpanded = attributes.get(SessionReader.SessionAttribute.EXPAND.getText());
        String colorScale = attributes.get(SessionReader.SessionAttribute.COLOR_SCALE.getText());
        String fontSizeString = attributes.get("fontSize");
        String showDataRangeString = attributes.get("showDataRange");

        if (colorScale != null) {
            ColorScale cs = ColorScaleFactory.getScaleFromString(colorScale);
            // This test should not be neccessary, refactor to eliminate it
            if (cs instanceof ContinuousColorScale) {
                this.setColorScale((ContinuousColorScale) cs);
            }
        }

        if (name != null && name.length() > 0) {
            setName(name);
        } else if (displayName != null && displayName.length() > 0) {
            setName(displayName);
        }

        // Set visibility
        if (isVisible != null) {
            if (isVisible.equalsIgnoreCase("true")) {
                setVisible(true);
            } else {
                setVisible(false);
            }
        }

        if (showDataRangeString != null) {
            try {
                showDataRange = Boolean.parseBoolean(showDataRangeString);
            } catch (Exception e) {
                log.error("Error parsing data range: " + showDataRangeString);
            }
        }

        // Set height
        if (height != null) {
            try {
                setHeight(Integer.parseInt(height));
            } catch (NumberFormatException e) {
                log.error("Error restoring track height: " + height);
            }
        }

        if (fontSizeString != null) {
            try {
                setFontSize(Integer.parseInt(fontSizeString));
            } catch (NumberFormatException e) {
                log.error("Error restoring font size: " + fontSizeString);
            }
        }

        // Set color
        if (color != null) {
            try {
                String[] rgb = color.split(",");
                int red = Integer.parseInt(rgb[0]);
                int green = Integer.parseInt(rgb[1]);
                int blue = Integer.parseInt(rgb[2]);
                setColor(new Color(red, green, blue));
            } catch (NumberFormatException e) {
                log.error("Error restoring color: " + color);
            }
        }

        // Set rendererClass
        if (rendererType != null) {
            Class rendererClass = RendererFactory.getRendererClass(rendererType);
            if (rendererClass != null) {
                setRendererClass(rendererClass);
            }
        }

        // Set window function
        if (windowFunction != null) {
            setStatType(WindowFunction.getWindowFunction(windowFunction));
        }

        // Set DataRange
        if (scale != null) {
            String[] axis = scale.split(",");
            float minimum = Float.parseFloat(axis[0]);
            float baseline = Float.parseFloat(axis[1]);
            float maximum = Float.parseFloat(axis[2]);
            setDataRange(new DataRange(minimum, baseline, maximum));
        }

        // Set expand
        if (isExpanded != null) {
            setExpanded(isExpanded.equalsIgnoreCase("true"));
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

    public String getId_142() {
        return id_142 == null ? getName() : id_142;
    }

    public void setId_142(String id_142) {
        this.id_142 = id_142;
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
            return (float) (Math.log(Math.max(Float.MIN_VALUE, value) / norm) / DataSourceTrack.log2);
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
}
