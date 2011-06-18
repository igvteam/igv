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
