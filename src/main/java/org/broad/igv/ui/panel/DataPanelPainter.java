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
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.sam.InsertionManager;
import org.broad.igv.sam.InsertionMarker;
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


        //

        Graphics2D graphics2D = context.getGraphics2D("BACKGROUND");
        graphics2D.setBackground(background);
        graphics2D.clearRect(visibleRect.x, visibleRect.y, visibleRect.width, visibleRect.height);
        graphics2D.setColor(Color.BLACK);


        final ReferenceFrame referenceFrame = context.getReferenceFrame();
        context.setInsertionMarkers(InsertionManager.getInstance().getInsertions(referenceFrame.getChrName(), referenceFrame.getOrigin(), referenceFrame.getEnd()));

        InsertionMarker i = InsertionManager.getInstance().getSelectedInsertion(referenceFrame.getChrName());


        if (i != null) {
            final double start = referenceFrame.getOrigin();
            try {
                double scale = referenceFrame.getScale();
                int p0 = 0, w;

                w = (int) ((i.position - start) / scale);
                paintSection(groups, context, p0, visibleRect.y, w, visibleRect.height, 0);

                i.pixelPosition = p0 + w;
                p0 += w;

                context.multiframe = true;
                referenceFrame.origin = i.position;
                w = (int) Math.ceil(i.size / context.getScale());
                paintExpandedInsertion(i, groups, context, p0, visibleRect.y, w, visibleRect.height);

                int g = (int) (i.size / scale);
                p0 += g;
                context.getReferenceFrame().origin = i.position;
                w = width - p0;
                if (w > 0) {
                    paintSection(groups, context, p0, visibleRect.y, width - p0, visibleRect.height, g);
                }
            } finally {
                referenceFrame.origin = start;
            }
        } else {
            paintFrame(groups, context, width, visibleRect);
        }

    }


    private void paintSection(Collection<TrackGroup> groups, RenderContext context, int px, int py, int w, int h, int delta) {

        context.clearGraphicsCache();
        Graphics2D dG = context.getGraphics();
        Rectangle dRect = new Rectangle(0, py, w, h);
        dG.translate(delta, 0);
        dG.setClip(dRect);
        context.translateX = px;

        paintFrame(groups, context, w, dRect);

    }


    private void paintFrame(Collection<TrackGroup> groups, RenderContext dContext, int width, Rectangle dRect) {
        int trackX = 0;
        int trackY = 0;

        for (Iterator<TrackGroup> groupIter = groups.iterator(); groupIter.hasNext(); ) {
            TrackGroup group = groupIter.next();

            if (dRect != null && (trackY > dRect.y + dRect.height)) {
                break;
            }

            if (group.isVisible()) {
                if (groups.size() > 1) {
                    final Graphics2D greyGraphics = dContext.getGraphic2DForColor(UIConstants.LIGHT_GREY);
                    greyGraphics.fillRect(0, trackY + 1, width, UIConstants.groupGap - 1);
                    trackY += UIConstants.groupGap;
                }

                // Draw a line just above group.
                if (group.isDrawBorder()) {
                    Graphics2D graphics2D = dContext.getGraphic2DForColor(Color.black);
                    graphics2D.drawLine(0, trackY - 1, width, trackY - 1);
                }

                List<Track> trackList = group.getVisibleTracks();
                synchronized (trackList) {
                    for (Track track : trackList) {
                        if (track == null) continue;
                        int trackHeight = track.getHeight();
                        if (dRect != null) {
                            if (trackY > dRect.y + dRect.height) {
                                break;
                            } else if (trackY + trackHeight < dRect.y) {
                                if (track.isVisible()) {
                                    trackY += trackHeight;
                                }
                                continue;
                            }
                        }


                        if (track.isVisible()) {
                            Rectangle rect = new Rectangle(trackX, trackY, width, trackHeight);
                            draw(track, rect, dContext);
                            trackY += trackHeight;
                        }
                    }
                }

                // Draw a line just below group.
                if (group.isDrawBorder()) {
                    Graphics2D graphics2D = dContext.getGraphic2DForColor(Color.black);
                    graphics2D.drawLine(0, trackY, width, trackY);
                }
            }
        }
    }


    private void paintExpandedInsertion(InsertionMarker insertionMarker, Collection<TrackGroup> groups, RenderContext context, int px, int py, int w, int h) {

        context.clearGraphicsCache();
        Graphics2D dG = context.getGraphics();
        Rectangle dRect = new Rectangle(0, py, w, h);
        dG.translate(px, 0);
        dG.setClip(dRect);
        context.translateX = px;

        int trackY = 0;

        for (Iterator<TrackGroup> groupIter = groups.iterator(); groupIter.hasNext(); ) {
            TrackGroup group = groupIter.next();
            if (group.isVisible()) {
                List<Track> trackList = group.getVisibleTracks();
                synchronized (trackList) {

                    for (Track track : trackList) {
                        if (track == null) continue;
                        int trackHeight = track.getHeight();
                        if (dRect != null) {
                            if (trackY > dRect.y + dRect.height) {
                                break;
                            } else if (trackY + trackHeight < dRect.y) {
                                if (track.isVisible()) {
                                    trackY += trackHeight;
                                }
                                continue;
                            }
                        }

                        if (track instanceof AlignmentTrack && track.isVisible()) {
                            Rectangle rect = new Rectangle(dRect.x, trackY, dRect.width, trackHeight);
                            ((AlignmentTrack) track).renderExpandedInsertion(insertionMarker, context, rect);
                        }

                        if (track.isVisible()) {
                            trackY += trackHeight;
                        }
                    }
                }
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


}


