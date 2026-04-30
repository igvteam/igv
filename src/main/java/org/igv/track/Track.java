package org.igv.track;


import htsjdk.tribble.Feature;
import htsjdk.tribble.NamedFeature;
import org.igv.Globals;
import org.igv.renderer.ContinuousColorScale;
import org.igv.renderer.DataRange;
import org.igv.renderer.Renderer;
import org.igv.sample.SampleGroup;
import org.igv.ui.panel.MouseableRegion;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.ui.panel.TrackPanelScrollPane;
import org.igv.util.ResourceLocator;
import org.json.JSONObject;
import org.w3c.dom.Element;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * @author jrobinso
 */
public interface Track {


    enum DisplayMode {
        COLLAPSED, SQUISHED, EXPANDED, FULL
    }

    /**
     * Return the igv.js compatible track type.  This is used for session export.
     *
     * @return
     */
    default TrackType getType() {
        return TrackType.annotation;
    }

    /**
     * Return an integer that defines the order of the track relative to other tracks.  Tracks are sorted by this value
     * in ascending order.  The default is 0, so tracks with a positive order will be below tracks with a negative order.
     * Note: Data type long is used to be compatible with igv.js, which uses MAX_SAFE_INTEGER, decidely not
     * safe for Java int.
     *
     * @return
     */
    default long getOrder() {
        return 0;
    }

    /**
     * Set the order of the track relative to other tracks.
     *
     * @param order the order value
     */
    default void setOrder(long order) {
        // default implementation does nothing
    }

    void setViewport(TrackPanelScrollPane trackPanel);

    TrackPanelScrollPane getViewport();

    /**
     * Return an identifier for the track.   The identifier should be unique in the context of a session.  For
     * files that produce a single track (e.g. wig) the absolute path name for the underlying file can be
     * used.
     *
     * @return
     */
    String getId();

    /**
     * Supports the whole genome view.  Default to true, but some tracks and data sources do not provide
     * whole genome features or scores.
     *
     * @return
     */
    default boolean supportsWholeGenome() {
        return true;
    }

    /**
     * Return true if the track is ready to paint (has all required data loaded).
     *
     * @param frame
     * @return
     */
    default boolean isReadyToPaint(ReferenceFrame frame) {
        return true;
    }

    /**
     * Load required resources to paint the reference frame.
     *
     * @param frame
     */
    default void load(ReferenceFrame frame) {
        // do nothing
    }

    /**
     * Return true if a track can be filtered by sample annotation.
     *
     * @return
     */
    default boolean isFilterable() {
        return true;
    }


    /**
     * Render the portion of the track overlapping visitbleRect.  In many cases visibleRect cover the entire track,
     * but in the case of scrolling it may not.
     *
     * @param context the render context
     */
    void render(RenderContext context);

    /**
     * Render the name of the track.
     *
     * @param graphics
     * @param trackRectangle
     * @param visibleRect
     */
    void renderName(Graphics2D graphics, Rectangle trackRectangle, Rectangle visibleRect);

    void renderAttributes(Graphics2D graphics, Rectangle trackRectangle,
                          List<String> names, List<MouseableRegion> mouseRegions);

    void setName(String name);

    String getName();

    String getDisplayName();

    default String getTooltipText(int y) {
        return getName();
    }

    ResourceLocator getResourceLocator();

    /**
     * Return ALL ResourceLocators associated with this track.
     * For most tracks, this will just be a collection of size 1
     *
     * @return
     */
    Collection<ResourceLocator> getResourceLocators();

    void setAttributeValue(String key, String value);

    void removeAttribute(String key);

    String getAttributeValue(String attributeKey);

    // Interface needed for tracks that contain multiple sub-tracks (samples).
    default String getAttribute(String name, String key) {
        return getAttributeValue(key);
    }

    void setVisible(boolean isVisible);

    boolean isVisible();

    default boolean isNumeric() {
        return false;
    }

    void setDataType(DataType type);

    DataType getDataType();

    void setY(int top);

    int getY();

    void setColorScale(ContinuousColorScale colorScale);

    ContinuousColorScale getColorScale();

    /**
     * Set the viewport height of the track.   This is the visible height of the track, not necessarily the content height.
     *
     * @param height
     */
    void setHeight(int height);

    /**
     * Return the visible (viewport) height of the tracks.  For some tracks this can differ from the content height,
     * in which case scrollbars will be shown.
     *
     * @return
     */
    int getHeight();

    /**
     * Return the minimum height for this track. Tracks should not be resized below this value.
     *
     * @return the minimum height in pixels
     */
    default int getMinimumHeight() {
        return 20;
    }

    /**
     * Return the content height of the track.  For tracks with variable numbers of rows this is calculated.
     *
     * @return
     */
    default int getContentHeight() {
        return 25;
    }

    /**
     * Return the height of a single row within a track. Applicable to tracks with multiple rows, including
     * alignment tracks, variant tracks, and feature tracks in "expand" mode.
     *
     * @return
     */
    int getRowHeight();

    void setRowHeight(int rowHeight);

    /**
     * Shrink the track to its minimum useful height. Default implementation sets the track height to its
     * minimum height. Tracks with variable content height (e.g. multi-row tracks) should override this to
     * use the larger of the content height or minimum height. Only shrinks; never grows the current height.
     */
    default void minimizeHeight() {
        setHeight(Math.min(getHeight(), getMinimumHeight()));
    }

    /**
     * Return the number of data rows currently displayed. Used to compute the row height
     * needed to fit all content into the viewport without scrolling.
     */
    default int getNumRows() {
        return 1;
    }

    /** Height in pixels reserved at the top of the track (not available for row content). */
    default int getReservedHeight() {
        return 0;
    }

    /**
     * Manually specify the data range.
     * {@code autoScale} must be turned off elsewhere, if applicable
     *
     * @param axisDefinition
     */
    void setDataRange(DataRange axisDefinition);

    DataRange getDataRange();

    Color getColor();

    default Color getDefaultColor() {
        return Globals.isDarkMode() ? Color.cyan : Color.blue.brighter();
    }

    void setColor(Color color);

    Color getAltColor();

    void setAltColor(Color color);

    void setWindowFunction(WindowFunction type);

    WindowFunction getWindowFunction();

    void setRendererClass(Class rc);

    Renderer getRenderer();

    boolean isShowDataRange();

    String getValueStringAt(String chr, double position, int mouseX, int mouseY, ReferenceFrame frame);

    float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, String frameName);

    float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, String frameName, List<Track> tracks);

    void setFontSize(int h);

    int getFontSize();

    boolean handleDataClick(TrackClickEvent e);

    void handleNameClick(MouseEvent e);

    Collection<WindowFunction> getAvailableWindowFunctions();

    void setProperties(TrackProperties trackProperties);

    Feature getFeatureAtMousePosition(TrackClickEvent e);

    void setSampleId(String sampleId);

    float logScaleData(float dataY);

    boolean isRegionScoreType(RegionScoreType type);

    int getVisibilityWindow();

    void setVisibilityWindow(int i);

    boolean isItemRGB();

    boolean isUseScore();

    default boolean hasDisplayMode() {
        return false;
    }

    DisplayMode getDisplayMode();

    void setDisplayMode(DisplayMode mode);

    /**
     * Return a list of menu items for the track's popup menu.  Separators can be
     * included with {@code new JPopupMenu.Separator()}.
     *
     * @param te
     * @return list of menu items, or null if no custom menu items
     */
    default List<Component> getPopupMenuItems(final TrackClickEvent te) {
        return null;
    }

    boolean isDrawYLine();

    float getYLine();

    void unload();

    boolean getAutoScale();

    void setAutoScale(boolean autoScale);

    default boolean isShowFeatureNames() {
        return true;
    }

    default void setShowFeatureNames(boolean b) {
    }

    /**
     * Return the java property or attribute for the feature display name.  Default is "null", in which case the
     * feature "name" property will be used.
     *
     * @return
     */
    default String getLabelField() {
        return null;
    }

    default void setTrackLine(String trackLine) {
    }

    /**
     * Return true if the track can be searched for a feature by name.
     *
     * @return
     */
    default boolean isSearchable() {
        return false;
    }

    default List<NamedFeature> search(String token) {
        return null;
    }


    String getSample();

    void repaint();

    /**
     * Return the number of samples in this track.  For most tracks this will be 1.
     *
     * @return
     */
    default int sampleCount() {
        return 1;
    }

    default List<SampleGroup> getSampleGroups() {
        return Collections.EMPTY_LIST;
    }

    default int getSampleHeight() {
        return getHeight();
    }

    default int getSampleOffset() {
        return 0;
    }

    default void sortSamplesByValue(String chr, int start, int end, RegionScoreType type) {
        // no op, override in subclass if needed
    }

    void setFeatureInfoURL(String featureInfoURL);

    /**
     * Restore object state from an XML element
     */

    default void unmarshalXML(Element element, Integer version) {
    }

    default void marshalJSON(JSONObject trackJson) {

    }

    default void unmarshalJSON(JSONObject jsonObject) {
    }


}
