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

import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.stream.Collectors;

/**
 * This class added to preload data when using gene lists.   Its actually more general than that,  but its motivation
 * stems from the need to provide a wait cursor when loading gene lists.
 *
 * @author jrobinso
 * @date Mar 22, 2011
 */
public class Preloader {

    public static synchronized CompletableFuture preload() {


        Collection<Track> trackList = visibleTracks();
        final List<ReferenceFrame> frames = FrameManager.getFrames();
        List<CompletableFuture> futures = new ArrayList(trackList.size() * frames.size());
        for (ReferenceFrame frame : frames) {
            for (Track track : trackList) {
                if (track.isReadyToPaint(frame) == false) {
                    futures.add(CompletableFuture.runAsync(() -> {
                        System.out.println("Preload " + Thread.currentThread());
                        System.out.println("Preloading" + track);
                        track.load(frame);
                    }));
                }
            }
        }

        if (futures.size() > 0) {
            WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();
            CompletableFuture[] futureArray = futures.toArray(new CompletableFuture[futures.size()]);
            return CompletableFuture.allOf(futureArray).thenRun(() -> WaitCursorManager.removeWaitCursor(token));
        }
        else {
            return CompletableFuture.completedFuture(null);
        }
    }

    public  static List<Track> visibleTracks() {
        return IGV.getInstance().getAllTracks().stream().
                filter(Track::isVisible).
                collect(Collectors.toList());
    }


}