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

import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.tribble.Feature;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

/**
 * @author jrobinso
 */
public class GeneTrack implements Track {

    SequenceTrack sequenceTrack;
    FeatureTrack featureTrack;

    /**
     * Constructs ...
     *
     * @param featureTrack
     * @param sequenceTrack
     */
    public GeneTrack(FeatureTrack featureTrack, SequenceTrack sequenceTrack) {
        this.featureTrack = featureTrack;
        this.sequenceTrack = sequenceTrack;
    }

    /**
     * Method description
     *
     * @param context
     * @param rect
     */
    public void render(RenderContext context, Rectangle rect) {

        // Split rects
        int seqHeight = sequenceTrack.getHeight();
        if (seqHeight > 0) {
            Rectangle seqRect = new Rectangle(rect);
            seqRect.height = seqHeight;
            sequenceTrack.render(context, seqRect);
        }

        rect.y += seqHeight;
        rect.height -= seqHeight;
        featureTrack.render(context, rect);
    }

    /**
     * Method description
     *
     * @param context
     * @param rect
     */
    public void renderBorder(RenderContext context, Rectangle rect) {

        // Split rects
        int seqHeight = sequenceTrack.getHeight();
        if (seqHeight > 0) {
            Rectangle seqRect = new Rectangle(rect);
            seqRect.height = seqHeight;
            sequenceTrack.renderBorder(context, seqRect);
        }

        rect.y += seqHeight;
        rect.height -= seqHeight;
        featureTrack.renderBorder(context, rect);
    }

    /**
     * Method description
     *
     * @param context
     * @param rect
     */
    public void renderAxis(RenderContext context, Rectangle rect) {

        // Split rects
        int seqHeight = sequenceTrack.getHeight();
        if (seqHeight > 0) {
            Rectangle seqRect = new Rectangle(rect);
            seqRect.height = seqHeight;
            sequenceTrack.renderAxis(context, seqRect);
        }

        rect.y += seqHeight;
        rect.height -= seqHeight;
        featureTrack.renderAxis(context, rect);
    }

    /**
     * Method description
     *
     * @param context
     * @param rect
     */
    public void overlay(RenderContext context, Rectangle rect) {
        featureTrack.render(context, rect);
        sequenceTrack.render(context, rect);
    }

    public String getId() {
        return featureTrack.getId();
    }

    /**
     * Method description
     *
     * @return
     */
    public String getName() {
        return featureTrack.getName();
    }

    public String getSampleId() {
        return featureTrack.getSampleId();
    }

    /**
     * Method description
     *
     * @return
     */
    public String getDisplayName() {
        return featureTrack.getName();
    }


    /**
     * Method description
     *
     * @return
     */
    public ResourceLocator getResourceLocator() {
        return featureTrack.getResourceLocator();
    }

    /**
     * Method description
     *
     * @param key
     * @param value
     */
    public void setAttributeValue(String key, String value) {
        featureTrack.setAttributeValue(key, value);
    }

    /**
     * Method description
     *
     * @param attributeKey
     * @return
     */
    public String getAttributeValue(String attributeKey) {
        return featureTrack.getAttributeValue(attributeKey);
    }

    /**
     * Method description
     *
     * @param isVisible
     */
    public void setVisible(boolean isVisible) {
        featureTrack.setVisible(isVisible);
        sequenceTrack.setVisible(isVisible);
    }

    /**
     * Method description
     *
     * @return
     */
    public boolean isVisible() {
        return featureTrack.isVisible() || sequenceTrack.isVisible();
    }

    /**
     * Method description
     *
     * @param type
     */
    public void setTrackType(TrackType type) {
        featureTrack.setTrackType(type);
    }

    /**
     * Method description
     *
     * @return
     */
    public TrackType getTrackType() {
        return featureTrack.getTrackType();
    }

    /**
     * Method description
     *
     * @param preferredHeight
     */
    public void setHeight(int preferredHeight) {
        featureTrack.setHeight(preferredHeight);
    }

    /**
     * Method description
     *
     * @return
     */
    public int getHeight() {
        return featureTrack.getHeight() + sequenceTrack.getHeight();
    }

    public int getMinimumHeight() {
        return featureTrack.getMinimumHeight() + sequenceTrack.getMinimumHeight();
    }

    /**
     * Method description
     *
     * @param axisDefinition
     */
    public void setDataRange(DataRange axisDefinition) {
        featureTrack.setDataRange(axisDefinition);
    }

    /**
     * Method description
     *
     * @return
     */
    public DataRange getDataRange() {
        return featureTrack.getDataRange();
    }

    /**
     * Method description
     *
     * @return
     */
    public Color getColor() {
        return featureTrack.getColor();
    }

    /**
     * Method description
     *
     * @return
     */
    public Color getAltColor() {
        return featureTrack.getAltColor();
    }

    /**
     * Method description
     *
     * @return
     */
    public Color getMidColor() {
        return featureTrack.getMidColor();
    }

    /**
     * Method description
     *
     * @param color
     */
    public void setColor(Color color) {
        featureTrack.setColor(color);
    }

    /**
     * Method description
     *
     * @param color
     */
    public void setAltColor(Color color) {
        featureTrack.setAltColor(color);
    }

    /**
     * Method description
     *
     * @param color
     */
    public void setMidColor(Color color) {
        featureTrack.setMidColor(color);
    }

    /**
     * Method description
     *
     * @param type
     */
    public void setStatType(WindowFunction type) {
        featureTrack.setStatType(type);
    }

    /**
     * Method description
     *
     * @return
     */
    public WindowFunction getWindowFunction() {
        return featureTrack.getWindowFunction();
    }

    /**
     * Method description
     *
     * @param rc
     */
    public void setRendererClass(Class rc) {
        featureTrack.setRendererClass(rc);
    }

    /**
     * Method description
     *
     * @return
     */
    public Renderer getRenderer() {
        return featureTrack.getRenderer();
    }

    /**
     * Method description
     *
     * @param selected
     */
    public void setSelected(boolean selected) {
        featureTrack.setSelected(selected);
    }

    /**
     * Method description
     *
     * @return
     */
    public boolean isSelected() {
        return featureTrack.isSelected();
    }

    /**
     * Method description
     *
     * @return
     */
    public boolean isDraggable() {
        return featureTrack.isDraggable();
    }

    /**
     * Method description
     *
     * @return
     */
    public boolean isLogNormalized() {
        return featureTrack.isLogNormalized();
    }

    public boolean isShowDataRange() {
        return featureTrack.isShowDataRange();
    }

    /**
     * Method description
     *
     * @param chr
     * @param position
     * @param y
     * @return
     */
    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {

        // TODO -- pass to seqeucne or feature track depending on y
        return featureTrack.getValueStringAt(chr, position, y, frame);
    }

    /**
     * Method description
     *
     * @param chr
     * @param start
     * @param end
     * @param zoom
     * @param type
     * @param frame
     * @return
     */
    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, ReferenceFrame frame) {
        return featureTrack.getRegionScore(chr, start, end, zoom, type, frame);
    }

    /**
     * Method description
     *
     * @param timestamp
     */
    public void refreshData(long timestamp) {
        featureTrack.refreshData(timestamp);
    }

    public void setFontSize(int h) {
        featureTrack.setFontSize(h);
    }

    public int getFontSize() {
        return featureTrack.getFontSize();
    }

    /**
     * Method description
     *
     * @param filename
     */
    public void setSourceFile(String filename) {
        featureTrack.setSourceFile(filename);
    }

    /**
     * Method description
     *
     * @return
     */
    public String getSourceFile() {
        return featureTrack.getSourceFile();
    }

    /**
     * @param value
     */
    public void setExpanded(boolean value) {
        featureTrack.setExpanded(value);
    }

    /**
     * @return
     */
    public boolean isExpanded() {
        return featureTrack.isExpanded();
    }


    public boolean handleClick(TrackClickEvent e) {
        return featureTrack.handleClick(e);
    }

    /**
     * Method description
     *
     * @return
     */
    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return featureTrack.getAvailableWindowFunctions();
    }

    /**
     * Method description
     *
     * @param trackProperties
     */
    public void setTrackProperties(TrackProperties trackProperties) {
        featureTrack.setTrackProperties(trackProperties);
    }

    /**
     * Method description
     *
     * @param overlayVisible
     */
    public void setOverlayVisible(boolean overlayVisible) {
        featureTrack.setOverlayVisible(overlayVisible);
    }

    /**
     * Method description
     *
     * @param name
     */
    public void setName(String name) {
        featureTrack.setName(name);
    }

    public void setUrl(String url) {
        featureTrack.setUrl(url);
    }

    public Feature getFeatureAtMousePosition(TrackClickEvent e) {
        return featureTrack.getFeatureAtMousePosition(e);
    }

    public String getId_142() {
        return featureTrack.getId_142();
    }

    public void setSampleId(String sampleId) {
        featureTrack.setSampleId(sampleId);
    }

    public float logScaleData(float dataY) {
        return featureTrack.logScaleData(dataY);
    }

    public boolean isRegionScoreType(RegionScoreType type) {
        return featureTrack.isRegionScoreType(type);
    }

    public int getVisibilityWindow() {
        return featureTrack.getVisibilityWindow();
    }

    public void setVisibilityWindow(int i) {
        featureTrack.setVisibilityWindow(i);
    }

    public Map<String, String> getPersistentState() {
        return featureTrack.getPersistentState();
    }

    public void restorePersistentState(Map<String, String> attributes) {
        featureTrack.restorePersistentState(attributes);
    }

    public String getUrl() {
        return featureTrack.getUrl();
    }

    public void renderName(Graphics2D graphics, Rectangle rect, Rectangle visibleRect) {
        featureTrack.renderName(graphics, rect, visibleRect);
    }

    public int getPreferredHeight() {
        return getHeight();
    }

    public void setTop(int top) {
        featureTrack.setTop(top);
    }

    public int getTop() {
        return featureTrack.getTop();
    }

    public void setColorScale(ContinuousColorScale colorScale) {
        featureTrack.setColorScale(colorScale);
    }

    //public String getActualName() {
    //    return featureTrack.getActualName();
    //}

    public ContinuousColorScale getColorScale() {
        return featureTrack.getColorScale();
    }

    public Feature nextFeature(String chr, double center, boolean forward, ReferenceFrame frame) throws IOException {
        return featureTrack.nextFeature(chr, center, forward, frame);
   }
}
