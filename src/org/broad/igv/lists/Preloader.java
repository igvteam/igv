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

package org.broad.igv.lists;

import org.apache.batik.bridge.CursorManager;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.ui.panel.DataPanel;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.LongRunningTask;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.Collectors;

/**
 * This class added to preload data when using gene lists.   Its actually more general than that,  but its motivation
 * stems from the need to provide a wait cursor when loading gene lists.
 *
 * @author jrobinso
 * @date Mar 22, 2011
 */
public class Preloader {

    private static Logger log = Logger.getLogger(Preloader.class);

    private static final ExecutorService threadExecutor = Executors.newFixedThreadPool(5);

    public static List<Track> visibleTracks() {
        return IGV.getInstance().getAllTracks().stream().
                filter(Track::isVisible).
                collect(Collectors.toList());
    }


    public static synchronized void load(final DataPanel dataPanel) {

        ReferenceFrame frame = dataPanel.getFrame();
        Collection<Track> trackList = dataPanel.visibleTracks();
        List<CompletableFuture> futures = new ArrayList(trackList.size());
        boolean batchLoaded = false;
        for (Track track : trackList) {
            if (track.isReadyToPaint(frame) == false) {
                final Runnable runnable = () -> {
                    //     log.info("Loading " + track.getName() + " " + frame.getFormattedLocusString());
                    track.load(frame);

                    //     log.info("Loaded " + track.getName() + " " + frame.getFormattedLocusString());
                };

                if(Globals.isBatch()) {
                    runnable.run();
                    batchLoaded = true;
                } else {
                    futures.add(CompletableFuture.runAsync(runnable, threadExecutor));
                }
            }
        }

        if (futures.size() > 0 || batchLoaded) {
            final CompletableFuture[] futureArray = futures.toArray(new CompletableFuture[futures.size()]);
            WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();
            CompletableFuture.allOf(futureArray).thenRun(() -> {

                //log.info("Call repaint " + dataPanel.hashCode() + " " + dataPanel.allTracksLoaded());
                dataPanel.loadInProgress = false;
                WaitCursorManager.removeWaitCursor(token);
                dataPanel.repaint();

            });
        }
    }

    private static SequenceTrack findSequenceTrack(Collection<Track> trackList) {

        for(Track t : trackList) {
            if(t instanceof SequenceTrack) {
                return (SequenceTrack) t;
            }
        }
        return null;

    }

}