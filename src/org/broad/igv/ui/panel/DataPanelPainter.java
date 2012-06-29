/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.exome.ExomeBlock;
import org.broad.igv.feature.exome.ExomeReferenceFrame;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.color.ColorUtilities;

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

            if (frame.isExomeMode()) {

                ExomeReferenceFrame exomeFrame = (ExomeReferenceFrame) frame;

                int blockGap = exomeFrame.getBlockGap();

                Rectangle panelClip = visibleRect;

                List<ExomeBlock> blocks = exomeFrame.getBlocks();
                int idx = exomeFrame.getFirstBlockIdx();
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
                        tmpFrame.setOrigin(b.getGenomeStart(), false);


                        RenderContext tmpContext = new RenderContextImpl(null, exomeGraphics, tmpFrame, rect);
                        preloadTracks(groups, tmpContext, rect);

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

    private void  preloadTracks(final Collection<TrackGroup> groups,
                               final RenderContext context,
                               final Rectangle visibleRect) {

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

        Runnable runnable = new Runnable() {
            public void run() {
                for (Track track : visibleTracks) {
                    RenderContextImpl newContext = new RenderContextImpl(null, null, context.getReferenceFrame(), null);
                    track.preload(newContext);
                }
            }
        };

        runnable.run();

//        Future future = LongRunningTask.submit(runnable);
//        try {
//            future.get();
//        } catch (InterruptedException e) {
//            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//        } catch (ExecutionException e) {
//            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//        }

//        Thread workerThread = new Thread(runnable);
//        workerThread.start();
//        try {
//            workerThread.join();
//        } catch (InterruptedException e) {
//            log.error("Preload thread was interrupted", e);
//        }
    }
}


