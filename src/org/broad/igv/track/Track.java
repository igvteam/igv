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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.track;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.tribble.Feature;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.session.Persistable;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.Collection;

/**
 * @author jrobinso
 */
public interface Track extends Persistable {

    enum DisplayMode {
        DENSE, SQUISH, PACK
    }

    ;

    /**
     * Render the track contents in the rectangle
     *
     * @param context
     * @param rect
     */
    public void render(RenderContext context, Rectangle rect);

    /**
     * Render track borders in the rectangle
     *
     * @param context
     * @param rect
     */
    public void renderBorder(RenderContext context, Rectangle rect);

    /**
     * Render the track border
     *
     * @param context
     * @param rect
     */
    public void renderAxis(RenderContext context, Rectangle rect);

    /**
     * Render the track name
     */

    public void renderName(Graphics2D graphics, Rectangle trackRectangle, Rectangle visibleRectangle);

    /**
     * Render this track as an overlay in the track rectangle
     *
     * @param context
     * @param rect
     */
    public void overlay(RenderContext context, Rectangle rect);

    /**
     * @return The unique identifier for the track.  This should only be used in saving/restoring sessions.
     */
    public String getId();

    /**
     * Get the display name of the track.
     *
     * @return
     */
    public String getName();

    /**
     * @return ID for use with the sample info table.  Not
     */
    public String getSampleId();


    public void setName(String name);

    public void setUrl(String url);

    public ResourceLocator getResourceLocator();

    public void setAttributeValue(String key, String value);

    public String getAttributeValue(String attributeKey);

    public void setVisible(boolean isVisible);

    public boolean isVisible();

    public void setOverlayVisible(boolean overlayVisible);

    public void setTrackType(TrackType type);

    public TrackType getTrackType();

    public void setHeight(int preferredHeight);

    public void setTop(int top);

    public int getTop();

    public void setColorScale(ContinuousColorScale colorScale);

    public ContinuousColorScale getColorScale();

    public int getHeight();

    public int getPreferredHeight();

    public int getMinimumHeight();

    public void setDataRange(DataRange axisDefinition);

    public DataRange getDataRange();

    public Color getColor();

    public void setColor(Color color);

    public Color getAltColor();

    public void setAltColor(Color color);

    public Color getMidColor();

    public void setMidColor(Color color);

    public void setStatType(WindowFunction type);

    public WindowFunction getWindowFunction();

    public void setRendererClass(Class rc);

    public Renderer getRenderer();

    public void setSelected(boolean selected);

    public boolean isSelected();

    public boolean isDraggable();

    public boolean isLogNormalized();

    public boolean isShowDataRange();

    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame);

    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, ReferenceFrame frame);

    public void refreshData(long timestamp);

    public void setFontSize(int h);

    public int getFontSize();

    /**
     * Record the file from which this track was created.
     *
     * @param filename
     */
    public void setSourceFile(String filename);

    /**
     * Return the filename from which this track was created.
     *
     * @return
     */
    public String getSourceFile();


    public void setExpanded(boolean value);

    public boolean isExpanded();

    public boolean handleDataClick(TrackClickEvent e);

    public void handleNameClick(MouseEvent e);

    public Collection<WindowFunction> getAvailableWindowFunctions();

    public void setTrackProperties(TrackProperties trackProperties);

    Feature getFeatureAtMousePosition(TrackClickEvent e);

    public String getId_142();

    void setSampleId(String sampleId);

    float logScaleData(float dataY);

    boolean isRegionScoreType(RegionScoreType type);

    int getVisibilityWindow();

    void setVisibilityWindow(int i);
}
