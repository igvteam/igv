package org.igv.ui.panel;

import org.igv.track.Track;
import org.igv.ui.IGV;

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
