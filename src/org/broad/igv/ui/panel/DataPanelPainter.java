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


package org.broad.igv.ui.panel;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;

import java.awt.*;
import java.util.*;
import java.util.List;


public class DataPanelPainter {

    private static Logger log = Logger.getLogger(DataPanelPainter.class);

    public synchronized void paint(Collection<TrackGroup> groups,
                                   RenderContext context,
                                   int width,
                                   Color background,
                                   Rectangle visibleRect) {

        Graphics2D graphics2D = null;

        try {
            graphics2D = (Graphics2D) context.getGraphics().create();
            graphics2D.setBackground(background);
            graphics2D.clearRect(visibleRect.x, visibleRect.y, visibleRect.width, visibleRect.height);
            graphics2D.setColor(Color.BLACK);

            paintFrame(groups, context, width, visibleRect);


        } finally {
            graphics2D.dispose();
        }
    }

    private void paintFrame(Collection<TrackGroup> groups,
                            RenderContext context,
                            int width,
                            Rectangle visibleRect) {


        int trackX = 0;
        int trackY = 0;

        groupAutoscale(groups, context);


        for (Iterator<TrackGroup> groupIter = groups.iterator(); groupIter.hasNext(); ) {
            TrackGroup group = groupIter.next();

            if (visibleRect != null && (trackY > visibleRect.y + visibleRect.height)) {
                break;
            }

            if (group.isVisible()) {
                if (groups.size() > 1) {
                    final Graphics2D greyGraphics = context.getGraphic2DForColor(UIConstants.LIGHT_GREY);
                    greyGraphics.fillRect(0, trackY + 1, width, UIConstants.groupGap - 1);
                    trackY += UIConstants.groupGap;
                }

                // Draw a line just above group.
                if (group.isDrawBorder()) {
                    Graphics2D graphics2D = context.getGraphic2DForColor(Color.black);
                    graphics2D.drawLine(0, trackY - 1, width, trackY - 1);
                }

                List<Track> trackList = group.getVisibleTracks();
                synchronized (trackList) {
                    for (Track track : trackList) {
                        if (track == null) continue;
                        int trackHeight = track.getHeight();
                        if (visibleRect != null) {
                            if (trackY > visibleRect.y + visibleRect.height) {
                                break;
                            } else if (trackY + trackHeight < visibleRect.y) {
                                if (track.isVisible()) {
                                    trackY += trackHeight;
                                }
                                continue;
                            }
                        }


                        if (track.isVisible()) {
                            Rectangle rect = new Rectangle(trackX, trackY, width, trackHeight);
                            draw(track, rect, context);
                            trackY += trackHeight;
                        }
                    }
                }

                // Draw a line just below group.
                if (group.isDrawBorder()) {
                    Graphics2D graphics2D = context.getGraphic2DForColor(Color.black);
                    graphics2D.drawLine(0, trackY, width, trackY);
                }
            }
        }
    }

    private void groupAutoscale(Collection<TrackGroup> groups, RenderContext context) {
        // TODO -- gene list mode!!!
        // Loop through all tracks extracting autoscale groups.
        // Get "in view scores" for each track in the group
        // Set data range for each track in the group and set  "autoScale" off (if its on).
        // Use  "groupAutoscale" tag?

        Map<String, List<Track>> autoscaleGroups = new HashMap<String, List<Track>>();
        for (Iterator<TrackGroup> groupIter = groups.iterator(); groupIter.hasNext(); ) {
            TrackGroup group = groupIter.next();
            List<Track> trackList = group.getVisibleTracks();
            synchronized (trackList) {
                for (Track track : trackList) {
                    String asGroup = track.getAttributeValue(AttributeManager.GROUP_AUTOSCALE);
                    if (asGroup != null) {
                        if (!autoscaleGroups.containsKey(asGroup)) {
                            autoscaleGroups.put(asGroup, new ArrayList<Track>());
                        }
                        autoscaleGroups.get(asGroup).add(track);
                    }
                }
            }
        }

        if (autoscaleGroups.size() > 0) {
            for (List<Track> tracks : autoscaleGroups.values()) {
                autoscale(context, tracks);
            }
        }
    }

    final private void draw(Track track, Rectangle rect, RenderContext context) {

        track.render(context, rect);

        // Get overlays

        List<Track> overlayTracks = IGV.getInstance().getOverlayTracks(track);
        if (overlayTracks != null) {
            for (Track overlayTrack : overlayTracks) {

                // Don't overlay on self
                if (overlayTrack != track) {
                    overlayTrack.overlay(context, rect);
                }
            }
        }

    }

    private List<Track> getVisibleTracks(final Collection<TrackGroup> groups) {
        // Find the tracks that need loaded, we go to this bother to avoid loading tracks scrolled out of view
        final List<Track> visibleTracks = new ArrayList<Track>();
        for (Iterator<TrackGroup> groupIter = groups.iterator(); groupIter.hasNext(); ) {
            TrackGroup group = groupIter.next();
            List<Track> trackList = new ArrayList(group.getVisibleTracks());
            for (Track track : trackList) {
                if (track != null && track.isVisible()) {
                    visibleTracks.add(track);
                }
            }
        }
        return visibleTracks;
    }


    private void autoscale(RenderContext context, List<Track> trackList) {
        int start = (int) context.getOrigin();
        int end = (int) context.getEndLocation() + 1;
        Rectangle visibleRect = context.getVisibleRect();
        List<LocusScore> inViewScores = new ArrayList<LocusScore>();
        synchronized (trackList) {
            for (Track track : trackList) {
                if (track instanceof DataTrack) {
                    inViewScores.addAll(((DataTrack) track).getInViewScores(context, visibleRect));
                }
            }

            if (inViewScores.size() > 0) {

                FeatureUtils.sortFeatureList(inViewScores);
                DataTrack.InViewInterval inter = DataTrack.computeScale(start, end, inViewScores);
                for (Track track : trackList) {
                    if (track instanceof DataTrack) {
                        DataRange dr = track.getDataRange();
                        float min = Math.min(0, inter.dataMin);
                        float base = Math.max(min, dr.getBaseline());
                        float max = inter.dataMax;
                        // Pathological case where min ~= max  (no data in view)
                        if (max - min <= (2 * Float.MIN_VALUE)) {
                            max = min + 1;
                        }

                        DataRange newDR = new DataRange(min, base, max, dr.isDrawBaseline());
                        newDR.setType(dr.getType());
                        track.setDataRange(newDR);
                    }
                }
            }
        }
    }


}


