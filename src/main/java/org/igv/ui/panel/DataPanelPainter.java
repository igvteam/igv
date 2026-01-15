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

    public static final Color GROUP_BORDER_COLOR = Globals.isDarkMode() ? Color.white : Color.black;
    public static final int TRACK_BORDER_HEIGHT = 3;
    private static Logger log = LogManager.getLogger(DataPanelPainter.class);

    public synchronized void paint(Collection<TrackGroup> groups,
                                   RenderContext context,
                                   Color background,
                                   Rectangle visibleRect) {

        Graphics2D borderGraphics = null;
        Color borderColor;
        List<Rectangle> borderRects = null;
        final boolean drawBorders = PreferencesManager.getPreferences().getAsBoolean(Constants.TRACK_DRAW_BORDERS);
        if (drawBorders) {
            borderGraphics = (Graphics2D) context.getGraphics().create();
            borderColor = PreferencesManager.getPreferences().getAsColor(Constants.TRACK_BORDER_COLOR);
            borderGraphics.setColor(borderColor);
        }

        try {
            Graphics2D graphics2D = context.getGraphics2D("BACKGROUND");
            graphics2D.setBackground(background);
            graphics2D.clearRect(visibleRect.x, visibleRect.y, visibleRect.width, visibleRect.height);

            final ReferenceFrame referenceFrame = context.getReferenceFrame();
            InsertionMarker i = referenceFrame.getExpandedInsertion(); //insertionManager.getSelectedInsertion(referenceFrame.getChrName());

            if (i != null) {

                Collection<AlignmentTrack> tracks = IGV.getInstance().getAlignmentTracks();
                int maxVizWindow = tracks.size() == 0 ? 0 : tracks.stream().mapToInt(t -> t.getVisibilityWindow()).max().getAsInt();

                final double start = referenceFrame.getOrigin();
                final double frameExtent = referenceFrame.getEnd() - referenceFrame.getOrigin();
                double scale = referenceFrame.getScale();
                int insertionPixelPosition = (int) ((i.position - start) / scale);
                int insertionPixelWidth = (int) Math.ceil(i.size / context.getScale());


                if (frameExtent < maxVizWindow &&
                        insertionPixelWidth > 1 &&
                        insertionPixelPosition < visibleRect.width &&
                        insertionPixelPosition + insertionPixelWidth > 0) {

                    context.expandedInsertionPosition = i.position;

                    // Paint section left of insertion
                    if (insertionPixelPosition > 0) {
                        int width = insertionPixelPosition;
                        RenderContext leftContext = shiftRenderContext(context, context.getOrigin(), 0, width);
                        borderRects = paintFrame(groups, leftContext);
                    }

                    // Paint expanded insertion
                    RenderContext insertionContext = shiftRenderContext(context, i.position, insertionPixelPosition, insertionPixelWidth);
                    paintExpandedInsertion(i, groups, insertionContext);

                    // Paint section to right of insertion
                    int p0 = insertionPixelPosition + insertionPixelWidth;
                    int w = visibleRect.width - p0;
                    if (w > 0) {
                        RenderContext rightContext = shiftRenderContext(context, i.position, p0, w);
                        rightContext.multiframe = true;
                        borderRects = paintFrame(groups, rightContext);
                    }

                } else {
                    // Insertion is out of view
                    borderRects = paintFrame(groups, context);
                }
            } else {
                // No expanded insertion
                borderRects = paintFrame(groups, context);
            }

            if (drawBorders) {
                for (Rectangle rectangle : borderRects) {
                    int bc = TRACK_BORDER_HEIGHT / 2;
                    borderGraphics.drawLine(0, rectangle.y + rectangle.height + bc, visibleRect.width, rectangle.y + rectangle.height + bc);
                }
            }

        } finally {
            if (context != null) {
                context.dispose();
            }
            if (borderGraphics != null) {
                borderGraphics.dispose();
            }
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


    private List<Rectangle> paintFrame(Collection<TrackGroup> groups, RenderContext dContext) {

        int trackX = 0;
        int trackY = 0;
        Rectangle dRect = dContext.visibleRect;
        List<Rectangle> borderRects = new ArrayList<>();

        for (Iterator<TrackGroup> groupIter = groups.iterator(); groupIter.hasNext(); ) {
            TrackGroup group = groupIter.next();

            if (dRect != null && (trackY > dRect.y + dRect.height)) {
                break;
            }

            if (group.isVisible()) {
                if (groups.size() > 1) {
                    final Graphics2D greyGraphics = dContext.getGraphic2DForColor(UIConstants.LIGHT_GREY);
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
                            borderRects.add(rect);
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
        return borderRects;
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


