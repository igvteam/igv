package org.igv.ui.panel;

import org.igv.Globals;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.sam.AlignmentTrack;
import org.igv.sam.InsertionMarker;
import org.igv.track.RenderContext;
import org.igv.track.Track;
import org.igv.track.TrackGroup;
import org.igv.ui.IGV;
import org.igv.ui.UIConstants;

import java.awt.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;


public class DataPanelPainter {

    public static final Color GROUP_BORDER_COLOR = Globals.isDarkMode() ? Color.darkGray : Color.lightGray;
    public static final int TRACK_BORDER_HEIGHT = 3;
    private static Logger log = LogManager.getLogger(DataPanelPainter.class);

    public synchronized void paint(Collection<TrackGroup> groups,
                                   RenderContext context,
                                   Color background,
                                   Rectangle visibleRect) {
        try {
            // Clear background
            Graphics2D graphics2D = context.getGraphics2D("BACKGROUND");
            graphics2D.setBackground(background);
            graphics2D.clearRect(visibleRect.x, visibleRect.y, visibleRect.width, visibleRect.height);

            // Paint tracks, handling expanded insertion if present
            InsertionMarker insertion = context.getReferenceFrame().getExpandedInsertion();
            if (shouldPaintExpandedInsertion(insertion, context, visibleRect)) {
                paintWithExpandedInsertion(groups, context, visibleRect, insertion);
            } else {
                paintFrame(groups, context);
            }

            // Draw track borders if enabled
            if (PreferencesManager.getPreferences().getAsBoolean(Constants.TRACK_DRAW_BORDERS)) {
                drawTrackBorders(groups, context);
            }
        } finally {
            context.dispose();
        }
    }

    private boolean shouldPaintExpandedInsertion(InsertionMarker insertion, RenderContext context, Rectangle visibleRect) {
        if (insertion == null) return false;

        Collection<AlignmentTrack> tracks = IGV.getInstance().getAlignmentTracks();
        if (tracks.isEmpty()) return false;

        ReferenceFrame frame = context.getReferenceFrame();
        int maxVizWindow = tracks.stream().mapToInt(Track::getVisibilityWindow).max().getAsInt();
        double frameExtent = frame.getEnd() - frame.getOrigin();
        if (frameExtent >= maxVizWindow) return false;

        double scale = frame.getScale();
        int pixelPos = (int) ((insertion.position - frame.getOrigin()) / scale);
        int pixelWidth = (int) Math.ceil(insertion.size / scale);

        return pixelWidth > 1 && pixelPos < visibleRect.width && pixelPos + pixelWidth > 0;
    }

    private void paintWithExpandedInsertion(Collection<TrackGroup> groups, RenderContext context,
                                            Rectangle visibleRect, InsertionMarker insertion) {
        ReferenceFrame frame = context.getReferenceFrame();
        double scale = frame.getScale();
        int pixelPos = (int) ((insertion.position - frame.getOrigin()) / scale);
        int pixelWidth = (int) Math.ceil(insertion.size / scale);

        context.expandedInsertionPosition = insertion.position;

        // Paint left of insertion
        if (pixelPos > 0) {
            paintFrame(groups, shiftRenderContext(context, context.getOrigin(), 0, pixelPos));
        }

        // Paint expanded insertion
        paintExpandedInsertion(insertion, groups, shiftRenderContext(context, insertion.position, pixelPos, pixelWidth));

        // Paint right of insertion
        int rightStart = pixelPos + pixelWidth;
        int rightWidth = visibleRect.width - rightStart;
        if (rightWidth > 0) {
            RenderContext rightContext = shiftRenderContext(context, insertion.position, rightStart, rightWidth);
            rightContext.multiframe = true;
            paintFrame(groups, rightContext);
        }
    }


    private RenderContext shiftRenderContext(RenderContext ctx, double position, int translateX, int pixelWidth) {
        RenderContext newContext = new RenderContext(ctx);
        newContext.getReferenceFrame().widthInPixels = pixelWidth;
        newContext.getReferenceFrame().origin = position;
        newContext.visibleRect = new Rectangle(0, ctx.visibleRect.y, pixelWidth, ctx.visibleRect.height);
        newContext.translateX = translateX;

        Graphics2D dG = newContext.getGraphics();
        dG.translate(translateX, 0);
        dG.setClip(newContext.visibleRect);

        return newContext;
    }


    private void paintFrame(Collection<TrackGroup> groups, RenderContext dContext) {

        int trackX = 0;
        int trackY = 0;
        Rectangle dRect = dContext.visibleRect;

        for (Iterator<TrackGroup> groupIter = groups.iterator(); groupIter.hasNext(); ) {
            TrackGroup group = groupIter.next();

            if (dRect != null && (trackY > dRect.y + dRect.height)) {
                break;
            }

            if (group.isVisible()) {
                if (groups.size() > 1) {
                    final Graphics2D greyGraphics = dContext.getGraphic2DForColor(GROUP_BORDER_COLOR);
                    greyGraphics.fillRect(0, trackY + 1, dRect.width, UIConstants.groupGap - 1);
                    trackY += UIConstants.groupGap;
                }

                // Draw a line just above group.
                if (group.isDrawBorder()) {
                    Graphics2D graphics2D = dContext.getGraphic2DForColor(GROUP_BORDER_COLOR);
                    graphics2D.drawLine(0, trackY - 1, dRect.width, trackY - 1);
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
                            Rectangle rect = new Rectangle(trackX, trackY, dRect.width, trackHeight);
                            Graphics2D g = dContext.getGraphics();
                            Shape originalClip = g.getClip();
                            // Clear cached graphics so new ones inherit the clip
                            dContext.clearGraphicsCache();
                            // Use clip() to intersect with existing clip, not setClip() which replaces it
                            g.clip(rect);
                            try {
                                draw(track, rect, dContext);
                            } catch (Exception e) {
                                log.error("Error rendering track: " + track.getName(), e);
                            } finally {
                                g.setClip(originalClip);
                                dContext.clearGraphicsCache();
                            }
                            trackY += trackHeight + TRACK_BORDER_HEIGHT; // Add 3 pixel gap between tracks
                        }
                    }
                }

                // Draw a line just below group.
                if (group.isDrawBorder()) {
                    Graphics2D graphics2D = dContext.getGraphic2DForColor(GROUP_BORDER_COLOR);
                    graphics2D.drawLine(0, trackY, dRect.width, trackY);
                }
            }
        }
    }


    private void drawTrackBorders(Collection<TrackGroup> groups, RenderContext context) {

        int trackY = 0;
        Rectangle dRect = context.getVisibleRect();
        int bottom = dRect.y + dRect.height;
        if (dRect == null) return;

        Graphics2D borderGraphics = (Graphics2D) context.getGraphics().create();
        Color borderColor = PreferencesManager.getPreferences().getAsColor(Constants.TRACK_BORDER_COLOR);
        borderGraphics.setColor(borderColor);

        try {
            for (Iterator<TrackGroup> groupIter = groups.iterator(); groupIter.hasNext(); ) {
                TrackGroup group = groupIter.next();
                if (group.isVisible()) {

                    List<Track> trackList = group.getVisibleTracks();

                    synchronized (trackList) {

                        for (Track track : trackList) {

                            if (track == null || !track.isVisible()) continue;
                            int trackHeight = track.getHeight();
                            if (trackY > bottom) {
                                // Beyond visible area, stop processing
                                break;
                            } else if (trackY + trackHeight < dRect.y) {
                                // Below visible area, just increment track y and continue
                                trackY += trackHeight + TRACK_BORDER_HEIGHT;
                            } else {
                                int borderY = trackY + trackHeight + TRACK_BORDER_HEIGHT / 2;
                                borderGraphics.drawLine(0, borderY, dRect.width, borderY);
                                trackY += trackHeight + TRACK_BORDER_HEIGHT;
                            }
                        }
                    }
                }
            }
        } finally {
            if (borderGraphics != null) {
                borderGraphics.dispose();
            }
        }
    }

    private void paintExpandedInsertion(InsertionMarker insertionMarker, Collection<TrackGroup> groups, RenderContext context) {


        int trackY = 0;
        Rectangle dRect = context.getVisibleRect();

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
                                // Below visible area, just increment track y and continue
                                if (track.isVisible()) {
                                    trackY += trackHeight + TRACK_BORDER_HEIGHT;
                                }
                                continue;
                            }
                        }

                        if (track instanceof AlignmentTrack && track.isVisible()) {
                            Rectangle rect = new Rectangle(dRect.x, trackY, dRect.width, trackHeight);
                            ((AlignmentTrack) track).renderExpandedInsertion(insertionMarker, context, rect);
                        }

                        if (track.isVisible()) {
                            trackY += trackHeight + TRACK_BORDER_HEIGHT;
                        }
                    }
                }
            }
        }
    }


    private void draw(Track track, Rectangle rect, RenderContext context) {

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


