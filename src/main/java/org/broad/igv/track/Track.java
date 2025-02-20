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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.track;


import htsjdk.tribble.Feature;
import htsjdk.tribble.NamedFeature;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.session.Persistable;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.MouseableRegion;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.Collection;
import java.util.List;

/**
 * @author jrobinso
 */
public interface Track extends Persistable {

    enum DisplayMode {
        COLLAPSED, SQUISHED, EXPANDED, FULL
    }

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
    boolean isReadyToPaint(ReferenceFrame frame);

    /**
     * Load required resources ot paint the reference frame.
     *
     * @param frame
     */
    void load(ReferenceFrame frame);

    /**
     * Return true if a track can be filtered by sample annotation.
     *
     * @return
     */
    default boolean isFilterable() {
        return true;
    }


    /**
     * Render the track in the supplied rectangle.  It is the responsibility of the track to draw within the
     * bounds of the rectangle.
     *
     * @param context the render context
     * @param rect    the track bounds, relative to the enclosing DataPanel bounds.
     */
    void render(RenderContext context, Rectangle rect);

    /**
     * Render the track as an overlay, presumably on another track.
     *
     * @param context the render context
     * @param rect    the track bounds, relative to the enclosing DataPanel bounds.
     */
    void overlay(RenderContext context, Rectangle rect);

    /**
     * Render the name of the track. Both the track and visible rectangles are supplied so the implementor
     * can adjust the placing of the name based on the current viewport.  This is used to center track names
     * on the viewport for large tracks that extend outside the viewport.
     *
     * @param graphics
     * @param trackRectangle   the track bounds, relative to the enclosing DataPanel bounds.
     * @param visibleRectangle
     */
    void renderName(Graphics2D graphics, Rectangle trackRectangle, Rectangle visibleRectangle);

    void renderAttributes(Graphics2D graphics, Rectangle trackRectangle, Rectangle visibleRect,
                          List<String> names, List<MouseableRegion> mouseRegions);

    void setName(String name);

    String getName();

    String getDisplayName();

    String getTooltipText(int y);

    String getSample();

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

    void setHeight(int preferredHeight);

    void setHeight(int preferredHeight, boolean force);

    void setY(int top);

    int getY();

    void setColorScale(ContinuousColorScale colorScale);

    ContinuousColorScale getColorScale();

    int getHeight();

    int getMinimumHeight();

    /**
     * Manually specify the data range.
     * {@code autoScale} must be turned off elsewhere, if applicable
     *
     * @param axisDefinition
     */
    void setDataRange(DataRange axisDefinition);

    DataRange getDataRange();

    Color getColor();

    default Color getExplicitColor() {
        return null;
    }

    void setColor(Color color);

    Color getAltColor();

    default Color getExplicitAltColor() {
        return null;
    }

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

    float getViewLimitMin();

    float getViewLimitMax();

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

    default void setTrackLine(String trackLine) {}

    /**
     * Return true if the track can be searched for a feature by name.
     *
     * @return
     */
    default boolean isSearchable() {
        return false;
    }

    default NamedFeature search(String token) {
        return null;
    }

    default void repaint() {
        IGV.getInstance().repaint(this);
    }

}
