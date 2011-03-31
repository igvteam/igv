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


import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.tribble.Feature;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.session.Persistable;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.Collection;

/**
 * @author jrobinso
 */
public interface Track extends Persistable {

    enum DisplayMode {
        COLLAPSED, SQUISHED, EXPANDED
    }

    // The unique identifier for the track.  This should only be used in saving/restoring sessions.
    public String getId();

    public void render(RenderContext context, Rectangle rect);

    public void renderName(Graphics2D graphics, Rectangle trackRectangle, Rectangle visibleRectangle);

    public void overlay(RenderContext context, Rectangle rect);

    void chromosomeChanged(String chrName);

    public void setName(String name);

    public String getName();

    public void setUrl(String url);

    public ResourceLocator getResourceLocator();

    public void setAttributeValue(String key, String value);

    public String getAttributeValue(String attributeKey);

    public void setVisible(boolean isVisible);

    public boolean isVisible();

    public void setOverlayVisible(boolean overlayVisible);

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

    public void setWindowFunction(WindowFunction type);

    public WindowFunction getWindowFunction();

    public void setRendererClass(Class rc);

    public Renderer getRenderer();

    public void setSelected(boolean selected);

    public boolean isSelected();

    public boolean isSortable();

    public boolean isLogNormalized();

    public boolean isShowDataRange();

    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame);

    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, ReferenceFrame frame);

    public void refreshData(long timestamp);

    public void setFontSize(int h);

    public int getFontSize();

    //public void setExpanded(boolean value);

    //public boolean isExpanded();

    public boolean handleDataClick(TrackClickEvent e);

    public void handleNameClick(MouseEvent e);

    public Collection<WindowFunction> getAvailableWindowFunctions();

    public void setTrackProperties(TrackProperties trackProperties);

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
 
    JPopupMenu getPopupMenu(final TrackClickEvent te);
    
}
