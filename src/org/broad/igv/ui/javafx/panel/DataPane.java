/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2017 Broad Institute
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
package org.broad.igv.ui.javafx.panel;

import org.broad.igv.track.Track;
import org.broad.igv.track.TrackGroup;
import org.broad.igv.ui.javafx.ResizableCanvas;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.util.List;
import java.util.stream.Collectors;

// Intended as the rough equivalent of the DataPanel class of the Swing UI.  Work in progress.
public class DataPane extends ResizableCanvas {

    private ReferenceFrame frame;
    private DataPaneContainer parent;

    public DataPane(ReferenceFrame frame, DataPaneContainer parent) {
        this.frame = frame;
        this.parent = parent;
    }

    public ReferenceFrame getFrame() {
        return frame;
    }

    // *** The following methods below this point copied over from TrackPanel as the functionality is the same. ***

    public boolean allTracksLoaded() {
        return parent.getTrackGroups().stream().
                filter(TrackGroup::isVisible).
                flatMap(trackGroup -> trackGroup.getVisibleTracks().stream()).
                allMatch(track -> track.isReadyToPaint(frame));
    }


    public List<Track> visibleTracks() {
        return parent.getTrackGroups().stream().
                filter(TrackGroup::isVisible).
                flatMap(trackGroup -> trackGroup.getVisibleTracks().stream()).
                collect(Collectors.toList());
    }

}
