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

package org.broad.igv.ui.dnd;

import org.broad.igv.track.Track;
import org.broad.igv.ui.panel.TrackPanel;

import java.util.ArrayList;
import java.util.List;
import java.awt.Point;


public class GhostDropEvent {
    private Point point;
    private Point startPoint;
    List<Track> tracks;
    List<TrackPanel> sourcePanels;
    boolean tracksDropped = false;


    public GhostDropEvent(Point startPoint, Point point, java.util.List<Track> tracks) {
        this.startPoint = startPoint;
        this.point = point;
        this.tracks = tracks;
        sourcePanels = new ArrayList();
    }

    public void addSourcePanel(TrackPanel panel) {
        sourcePanels.add(panel);
    }

    public void setTracksDropped(boolean bool) {
        tracksDropped = bool;
    }

    public boolean isTracksDropped() {
        return tracksDropped;
    }

    /**
     * Called when the destination panel is found.  If the tracks are dropped outside a valid destination panel
     * this is never called
     */
    public void removeTracksFromSource() {

        for(TrackPanel panel : sourcePanels) {
            panel.removeTracks(tracks);
        }
        sourcePanels.clear();

    }

    public Point getStartLocation() {
        return startPoint;
    }

    public Point getDropLocation() {
        return point;
    }

    public java.util.List<Track> getTracks() {
        return tracks;
    }
}
