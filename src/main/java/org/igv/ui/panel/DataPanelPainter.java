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
import org.igv.ui.IGV;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.util.Collection;


public class DataPanelPainter {

    public static final Color GROUP_BORDER_COLOR = Globals.isDarkMode() ? Color.darkGray : Color.lightGray;
    public static final int TRACK_BORDER_HEIGHT = 3;
    private static Logger log = LogManager.getLogger(DataPanelPainter.class);

    public synchronized void paint(Track track, RenderContext context) {
        try {

            // Paint tracks, handling expanded insertion if present
            InsertionMarker insertion = context.getReferenceFrame().getExpandedInsertion();
            if (shouldPaintExpandedInsertion(insertion, context)) {
                paintWithExpandedInsertion(track, context, insertion);
            } else {
                paintFrame(track, context);
            }

        } finally {
            context.dispose();
        }
    }

    private boolean shouldPaintExpandedInsertion(InsertionMarker insertion, RenderContext context) {
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

        return pixelWidth > 1 && pixelPos < context.getTrackRectangle().width && pixelPos + pixelWidth > 0;
    }

    private void paintWithExpandedInsertion(Track track, RenderContext context, InsertionMarker insertion) {

        ReferenceFrame frame = context.getReferenceFrame();
        double scale = frame.getScale();
        int pixelPos = (int) ((insertion.position - frame.getOrigin()) / scale);

        context.expandedInsertionPosition = insertion.position;

        // Paint left of insertion
        if (pixelPos > 0) {
            paintFrame(track, shiftRenderContext(context, context.getOrigin(), 0, pixelPos));
        }

        // Paint expanded insertion
        int pixelWidth = (int) Math.ceil(insertion.size / scale);
        paintExpandedInsertion(insertion, track, shiftRenderContext(context, insertion.position, pixelPos, pixelWidth));

        // Paint right of insertion
        int rightStart = pixelPos + pixelWidth;
        int rightWidth = context.getTrackRectangle().width - rightStart;
        if (rightWidth > 0) {
            RenderContext rightContext = shiftRenderContext(context, insertion.position, rightStart, rightWidth);
            rightContext.multiframe = true;
            paintFrame(track, rightContext);
        }
    }


    private RenderContext shiftRenderContext(RenderContext ctx, double position, int translateX, int pixelWidth) {
        RenderContext newContext = new RenderContext(ctx);
        newContext.getReferenceFrame().widthInPixels = pixelWidth;
        newContext.getReferenceFrame().origin = position;
        newContext.trackRectangle = new Rectangle(0, ctx.trackRectangle.y, pixelWidth, ctx.trackRectangle.height);
        newContext.translateX = translateX;

        Graphics2D dG = newContext.getGraphics();
        dG.translate(translateX, 0);
        dG.setClip(newContext.trackRectangle);

        return newContext;
    }


    private void paintFrame(Track track, RenderContext context) {

        if (track.isVisible()) {
            try {
//                float rowHeight = track.getRowHeight();
//                if (rowHeight > 0 && rowHeight != Math.round(rowHeight)) {
//                    paintScaled(track, context);
//                } else {
                    track.render(context);
 //               }
            } catch (Exception e) {
                log.error("Error rendering track: " + track.getName(), e);
            }
        }
    }

//    private void paintScaled(Track track, RenderContext context) {
//        Rectangle trackRect = context.getTrackRectangle();
//        int width = trackRect.width;
//        int imageHeight = track.getContentHeight();
//        if (imageHeight <= 0 || imageHeight == trackRect.height) {
//            track.render(context);
//            return;
//        }
//
//        Color bg = PreferencesManager.getPreferences().getAsColor(Constants.BACKGROUND_COLOR);
//        if (bg == null) bg = Color.WHITE;
//
//        BufferedImage image = new BufferedImage(width, imageHeight, BufferedImage.TYPE_INT_RGB);
//        Graphics2D ig = image.createGraphics();
//        ig.setColor(bg);
//        ig.fillRect(0, 0, width, imageHeight);
//
//        Rectangle offRect = new Rectangle(0, 0, width, imageHeight);
//        RenderContext offContext = new RenderContext(null, ig, context.getReferenceFrame(), offRect, offRect, offRect);
//        try {
//            track.render(offContext);
//        } finally {
//            offContext.dispose();
//            ig.dispose();
//        }
//
//        Graphics2D g = context.getGraphics();
//        Object prevHint = g.getRenderingHint(RenderingHints.KEY_INTERPOLATION);
//        g.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);
//        g.drawImage(image, trackRect.x, trackRect.y, trackRect.width, trackRect.height, null);
//        if (prevHint != null) {
//            g.setRenderingHint(RenderingHints.KEY_INTERPOLATION, prevHint);
//        }
//    }


    private void paintExpandedInsertion(InsertionMarker insertionMarker, Track track, RenderContext context) {
        if (track instanceof AlignmentTrack && track.isVisible()) {
            Rectangle dRect = context.getVisibleRect();

            ((AlignmentTrack) track).renderExpandedInsertion(insertionMarker, context, dRect);
        }

    }

}


