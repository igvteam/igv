/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.track;


import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.session.Persistable;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.MouseableRegion;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.Feature;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.Collection;
import java.util.List;

/**
 * @author jrobinso
 */
public interface Track extends Persistable{

    enum DisplayMode {
        COLLAPSED, SQUISHED, EXPANDED, ALTERNATIVE_SPLICE
    }

    /**
     * Return an identifier for the track.   The identifier should be unique in the context of a session.  For
     * files that produce a single track (e.g. wig) the absolute path name for the underlying file can be
     * used.
     *
     * @return
     */
    String getId();


    void load(RenderContext context);

    /**
     * Return true if a track can be filtered by sample annotation.
     *
     * @return
     */
    boolean isFilterable();


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
     * We store the y-coordinate at which this track was last rendered,
     * to avoid repeating borders/scales/etc. This resets that value
     */
    void resetLastY();

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

    String getNameValueString(int y);

    String getSample();

    void setUrl(String url);

    ResourceLocator getResourceLocator();

    /**
     * Return ALL ResourceLocators associated with this track.
     * For most tracks, this will just be a collection of size 1
     * @return
     */
    Collection<ResourceLocator> getResourceLocators();

    void setAttributeValue(String key, String value);

    String getAttributeValue(String attributeKey);

    void setVisible(boolean isVisible);

    boolean isVisible();

    void setOverlayed(boolean overlayVisible);

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
     * @param axisDefinition
     */
    void setDataRange(DataRange axisDefinition);

    boolean hasDataRange();

    DataRange getDataRange();

    Color getColor();

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

    String getValueStringAt(String chr, double position, int y, ReferenceFrame frame);

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

    void dispose();
    
    boolean getAutoScale();

    void setAutoScale(boolean autoScale);

}
