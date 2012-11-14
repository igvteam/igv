/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General  License (LGPL),
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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.track;


import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.MouseableRegion;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.tribble.Feature;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.session.Persistable;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 */
public interface Track extends Persistable {

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


    void preload(RenderContext context);

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

    String getNameValueString(int y);

    String getSample();

    void setUrl(String url);

    ResourceLocator getResourceLocator();

    void setAttributeValue(String key, String value);

    String getAttributeValue(String attributeKey);

    void setVisible(boolean isVisible);

    boolean isVisible();

    void setOverlayed(boolean overlayVisible);

    TrackType getTrackType();

    void setHeight(int preferredHeight);

    void setY(int top);

    int getY();

    void setColorScale(ContinuousColorScale colorScale);

    ContinuousColorScale getColorScale();

    int getHeight();

    int getMinimumHeight();

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

    void updateGenome(Genome genome);

}
