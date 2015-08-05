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

import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;

import java.awt.*;
import java.util.HashSet;
import java.util.Set;

/**
 * Regions that represent tracks
 */
public class MouseableRegion {

    private Shape region;
    private String text;
    private TrackCltn trackCltn;


    public MouseableRegion(Shape region, Track track) {

        this.region = region;

        if (track.getName().equals(track.getName())) {
            this.text = track.getName();
        } else {
            this.text = "<html>" + track.getName() + "<br>" + track.getName();
        }
        trackCltn = new SingleTrackRegion(track);
    }

    public MouseableRegion(Shape region, String name, String value) {

        this.region = region;
        this.text = (name.length() > 0 ? name + " = " : "") + value;
        trackCltn = new AttributePanelRegion(name, value);
    }

    public boolean containsPoint(double x, double y) {
        return region.contains(x, y);
    }


    public Rectangle getBounds() {
        return region.getBounds();
    }


    public String getText() {
        return text;
    }

    public Set<Track> getTracks() {
        return trackCltn.getTracks();
    }


    @Override
    public String toString() {
        return region.toString();
    }

    interface TrackCltn {

        public Set<Track> getTracks();

        public boolean contains(Track t);
    }

    class SingleTrackRegion implements TrackCltn {

        private Set<Track> tracks;

        public SingleTrackRegion(Track track) {
            tracks = new HashSet<Track>();
            tracks.add(track);
        }

        public Set<Track> getTracks() {
            return tracks;
        }

        public boolean contains(Track t) {
            return tracks.contains(t);
        }
    }

    class AttributePanelRegion implements TrackCltn {

        private String key;
        private String value;


        public AttributePanelRegion(String key, String value) {
            this.key = key.toUpperCase();
            this.value = value;
        }

        public Set<Track> getTracks() {
            Set<Track> selectedTracks = new HashSet();
            for (Track track : IGV.getInstance().getAllTracks()) {
                String attributeValue = track.getAttributeValue(key);
                if (attributeValue == null) {
                    continue;
                }
                if (attributeValue.equals(value)) {
                    selectedTracks.add(track);
                }
            }
            return selectedTracks;
        }

        public boolean contains(Track t) {
            return getTracks().contains(t);
        }
    }
}
