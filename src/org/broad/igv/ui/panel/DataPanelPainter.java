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
import org.broad.igv.feature.exome.ExomeBlock;
import org.broad.igv.feature.exome.ExomeReferenceFrame;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.RenderContextImpl;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackGroup;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;

import java.awt.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;


public class DataPanelPainter {

    private static Logger log = Logger.getLogger(DataPanelPainter.class);

    private static Color exomeBorderColor = new Color(190, 190, 255);

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

            ReferenceFrame frame = context.getReferenceFrame();

            resetLastY(groups);

            if (frame.isExomeMode()) {

                ExomeReferenceFrame exomeFrame = (ExomeReferenceFrame) frame;
                int blockGap = exomeFrame.getBlockGap();

                Rectangle panelClip = visibleRect;

                List<ExomeBlock> blocks = exomeFrame.getBlocks();
                int idx = exomeFrame.getFirstBlockIdx();

                RenderContext exomeContext = new RenderContextImpl(null, null, frame, visibleRect);
                preloadTracks(groups, exomeContext, visibleRect);

                ExomeBlock b;
                int lastPStart = 0;
                int pStart;
                int pEnd;
                int exomeOrigin = ((ExomeReferenceFrame) frame).getExomeOrigin();
                int visibleBlockCount = 0;

                do {
                    b = blocks.get(idx);

                    pStart = (int) ((b.getExomeStart() - exomeOrigin) / frame.getScale()) + visibleBlockCount * blockGap;
                    pEnd = (int) ((b.getExomeEnd() - exomeOrigin) / frame.getScale()) + visibleBlockCount * blockGap;

                    if (pEnd > lastPStart) {

                        lastPStart = pStart;
                        // Don't draw over previously drawn region -- can happen when zoomed out.


                        if (pEnd == pStart) pEnd++;


                        b.setScreenBounds(pStart, pEnd);

                        Rectangle rect = new Rectangle(pStart, visibleRect.y, pEnd - pStart, visibleRect.height);

                        Graphics2D exomeGraphics = (Graphics2D) context.getGraphics().create();

                        //Shape clip = exomeGraphics.getClip();

//                         Color c = ColorUtilities.randomColor(idx);
//                         exomeGraphics.setColor(c);
//                         exomeGraphics.fill(rect);
//                         exomeGraphics.setColor(Color.black);
//                         GraphicUtils.drawCenteredText(String.valueOf(idx), rect, exomeGraphics);

                        exomeGraphics.setClip(rect.intersection(panelClip));
                        exomeGraphics.translate(pStart, 0);

                        ReferenceFrame tmpFrame = new ReferenceFrame(frame);
                        tmpFrame.setOrigin(b.getGenomeStart());


                        RenderContext tmpContext = new RenderContextImpl(null, exomeGraphics, tmpFrame, rect);
                        paintFrame(groups, tmpContext, rect.width, rect);

                        tmpContext.dispose();
                        exomeGraphics.dispose();
                        visibleBlockCount++;
                    }
                    idx++;

                }
                while ((pStart < visibleRect.x + visibleRect.width) && idx < blocks.size());

                // Draw lines @ gene boundaries
                String chr = frame.getChrName();
                List<ExomeReferenceFrame.Gene> genes = exomeFrame.getGenes(chr);

                idx = FeatureUtils.getIndexBefore(frame.getOrigin(), genes);

                exomeOrigin = ((ExomeReferenceFrame) frame).getExomeOrigin();
                int top = visibleRect.y;

                int lastXDrawn = -1;
                Graphics2D lineGraphics = context.getGraphic2DForColor(exomeBorderColor);
                do {
                    ExomeReferenceFrame.Gene gene = genes.get(idx);
                    double exomeStart = exomeFrame.genomeToExomePosition(gene.getStart());
                    double exomeEnd = exomeFrame.genomeToExomePosition(gene.getEnd());
                    pStart = (int) ((exomeStart - exomeOrigin) / frame.getScale()) + visibleBlockCount * blockGap;
                    pEnd = (int) ((exomeEnd - exomeOrigin) / frame.getScale()) + visibleBlockCount * blockGap;

                    if (pStart != lastXDrawn) {
                        lineGraphics.drawLine(pStart, top, pStart, top + visibleRect.height);
                    }
                    lineGraphics.drawLine(pEnd, top, pEnd, top + visibleRect.height);
                    lastXDrawn = pEnd;
                    idx++;

                }
                while ((pStart < visibleRect.x + visibleRect.width) && idx < genes.size());


            } else {
                paintFrame(groups, context, width, visibleRect);
            }


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

                List<Track> trackList = group.getTracks();
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

    private void resetLastY(Collection<TrackGroup> groups) {
        for(TrackGroup group: groups){
            for(Track track: group.getTracks()){
                track.resetLastY();
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
            List<Track> trackList = new ArrayList(group.getTracks());
            for (Track track : trackList) {
                if (track != null && track.isVisible()) {
                    visibleTracks.add(track);
                }
            }
        }
        return visibleTracks;
    }

    private void preloadTracks(final Collection<TrackGroup> groups,
                               final RenderContext context,
                               final Rectangle visibleRect) {


        final List<Track> visibleTracks = getVisibleTracks(groups);

        for (Track track : visibleTracks) {
            RenderContextImpl newContext = new RenderContextImpl(null, null, context.getReferenceFrame(), visibleRect);
            track.load(newContext);
        }

    }
}


