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


package org.broad.igv.ui.panel;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackGroup;
import org.broad.igv.ui.IGVMainFrame;
import org.broad.igv.ui.UIConstants;

import java.awt.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

/**
 * @author jrobinso
 */
public class DataPanelPainter {

    private static Logger log = Logger.getLogger(DataPanelPainter.class);

    public void paint(Collection<TrackGroup> groups,
                      RenderContext context,
                      int width,
                      int height,
                      Color background,
                      Rectangle visibleRect) {

        Graphics2D graphics2D = null;

        try {
            graphics2D = (Graphics2D) context.getGraphics().create();

            graphics2D.setBackground(background);
            graphics2D.clearRect(0, 0, width, height);
            graphics2D.setColor(Color.BLACK);

            final Graphics2D greyGraphics = context.getGraphic2DForColor(UIConstants.VERY_LIGHT_GRAY);
            int trackX = 0;
            int trackY = 0;


            for (Iterator<TrackGroup> groupIter = groups.iterator(); groupIter.hasNext();) {
                TrackGroup group = groupIter.next();

                if (visibleRect != null && (trackY > visibleRect.y + visibleRect.height)) {
                    break;
                }

                if (group.isVisible()) {
                    if (groups.size() > 1) {
                        greyGraphics.fillRect(0, trackY + 1, width, UIConstants.groupGap - 1);
                        trackY += UIConstants.groupGap;
                    }

                    // Draw a line just above group.
                    if (group.isDrawBorder()) {
                        graphics2D.drawLine(0, trackY - 1, width, trackY - 1);
                    }

                    List<Track> trackList = group.getTracks();
                    for (Track track : trackList) {
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
                            trackY = draw(track, trackX, trackY, width, trackHeight, context);
                        }
                    }

                    // Draw a line just below group.
                    if (group.isDrawBorder()) {
                        graphics2D.drawLine(0, trackY, width, trackY);
                    }
                }
            }

            // Now draw borders.  This is done in a second loop to prevent data from
            // one track from overwriting border from an adjacent track
            trackX = 0;
            trackY = 0;
            for (Iterator<TrackGroup> groupIter = groups.iterator(); groupIter.hasNext();) {

                TrackGroup group = groupIter.next();
                if (group.isVisible()) {
                    if (groups.size() > 1) {
                        trackY += UIConstants.groupGap;
                    }

                    // Draw borders  -- copy track list to prevent concurrent modification exception
                    // Copy tracks to prevent concurrent modfication exception
                    List<Track> trackList = new ArrayList(group.getTracks());
                    for (Track track : trackList) {
                        if (track.isVisible()) {
                            Rectangle rect = new Rectangle(trackX, trackY, width, track.getHeight());
                            track.renderBorder(context, rect);
                            track.renderAxis(context, rect);
                            trackY += track.getHeight();
                        }
                    }
                }
            }
        }
        finally {
            graphics2D.dispose();
        }
    }

    final private int draw(Track track, int trackX, int trackY, int trackWidth, int trackHeight, RenderContext context) {


        // The rectangle in which the track is drawn.  The height is reduced
        Rectangle rect = new Rectangle(trackX, trackY, trackWidth, trackHeight);


        track.render(context, rect);

        // Get overlays

        List<Track> overlayTracks = IGVMainFrame.getInstance().getTrackManager().getOverlayTracks(track);
        if (overlayTracks != null) {
            for (Track overlayTrack : overlayTracks) {

                // Don't overlay on self
                if (overlayTrack != track) {
                    overlayTrack.overlay(context, rect);
                }
            }
        }

        // Change Y to point to the start of the next track
        trackY += trackHeight;

        return trackY;

    }
}
