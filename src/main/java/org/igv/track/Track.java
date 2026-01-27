package org.igv.track;


import htsjdk.tribble.Feature;
import htsjdk.tribble.NamedFeature;
import org.igv.Globals;
import org.igv.renderer.ContinuousColorScale;
import org.igv.renderer.DataRange;
import org.igv.renderer.Renderer;
import org.igv.session.Persistable;
import org.igv.ui.IGV;
import org.igv.ui.panel.*;
import org.igv.util.ResourceLocator;
import org.igv.util.TrackFilter;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;

/**
 * @author jrobinso
 */
public interface Track extends Persistable, AttributeSupplier {

    enum DisplayMode {
        COLLAPSED, SQUISHED, EXPANDED, FULL
    }

    void setViewport(TrackPanelScrollPane trackPanel);

    /**
     * Return an identifier for the track.   The identifier should be unique in the context of a session.  For
     * files that produce a single track (e.g. wig) the absolute path name for the underlying file can be
     * used.
     *
     * @return
     */
    String getId();

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
     * @param visibleRect
     */
    void render(RenderContext context, Rectangle visibleRect);

    /**
     * Render the track as an overlay, presumably on another track.
     *
     * @param context the render context
     * @param rect    the track bounds, relative to the enclosing DataPanel bounds.
     */
    default void overlay(RenderContext context, Rectangle rect) {
        // do nothing, overlay is optional
    }

    /**
     * Render the name of the track.
     *
     * @param graphics
     * @param trackRectangle
     */
    void renderName(Graphics2D graphics, Rectangle trackRectangle);

    void renderAttributes(Graphics2D graphics, Rectangle trackRectangle, Rectangle visibleRect,
                          List<String> names, List<MouseableRegion> mouseRegions);

    void setName(String name);

    String getName();

    String getDisplayName();

    default String getTooltipText(int y) {
        return getName();
    }

    String getSample();

    /**
     * Return the number of samples in this track.  For most tracks this will be 1.
     *
     * @return
     */
    default int sampleCount() {
        return 1;
    }

    default void sortSamplesByAttribute(Comparator<String> comparator) {
        // no op, override in subclass if needed
    }

    default void sortSamplesByValue(String chr, int start, int end, RegionScoreType type) {
        // no op, override in subclass if needed
    }

    default void filterSamples(TrackFilter trackFilter) {
        // no op, override in subclass if needed
    }

    void setFeatureInfoURL(String featureInfoURL);

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

    default boolean isAlignment() {
        return false;
    }

    void setOverlayed(boolean overlayVisible);

    void setTrackType(TrackType type);

    TrackType getTrackType();


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
     * Return the content height of the track.  For tracks with variable numbers of rows this is calculated.
     *
     * @return
     */
    default int getContentHeight() {
        return 20;
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

    void setSelected(boolean selected);

    boolean isSelected();

    boolean isSortable();

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

    DisplayMode getDisplayMode();

    void setDisplayMode(DisplayMode mode);

    IGVPopupMenu getPopupMenu(final TrackClickEvent te);

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

    default void groupSamplesByAttribute(String attributeKey) {
        // no op, override in subclass if needed
    }

    default void repaint() {
        if(IGV.hasInstance()) {
            IGV.getInstance().repaint(this);
        }
    }

}
