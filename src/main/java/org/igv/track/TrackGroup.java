package org.igv.track;

//~--- JDK imports ------------------------------------------------------------

import org.igv.logging.*;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.renderer.GraphicUtils;
import org.igv.ui.FontManager;
import org.igv.ui.IGV;

import java.awt.*;
import java.util.*;
import java.util.List;

/**
 * Container for a group of tracks.  Behaves as a single unit when sorting
 * by region score.
 *
 * @author jrobinso
 */
public class TrackGroup {

    private static Logger log = LogManager.getLogger(TrackGroup.class);

    /**
     * Key used to group tracks (e.g. SAMPLE_ID).
     */
    private String name;

    private List<Track> tracks;

    private boolean autoScale = false;

    public TrackGroup(String name) {
        this.name = name;
        tracks = Collections.synchronizedList(new ArrayList<Track>());
    }

    public boolean contains(Track track) {
        return tracks.contains(track);
    }

    public List<Track> getTracks() {
        return tracks;
    }

    public List<Track> getVisibleTracks() {
        List<Track> visibleTracks = new ArrayList<Track>();
        for(Track t : tracks) {
            if(t.isVisible()) visibleTracks.add(t);
        }
        return visibleTracks;
    }

    public boolean isAutoScale() {
        return autoScale;
    }

    public void setAutoScale(boolean autoScale) {
        this.autoScale = autoScale;
    }

    public int indexOf(Track track) {
        return tracks.indexOf(track);
    }


    public int size() {
        return tracks.size();
    }


    public void add(Track track) {
        if (track == null) {
            log.warn("Attempt to add null track");
        } else {
            tracks.add(track);
        }
    }

    public void add(int pos, Track track) {
        if (track == null) {
            log.warn("Attempt to add null track");
        } else {
            tracks.add(pos, track);
        }
    }

    public void addAll(Collection<Track> trackList) {
        tracks.addAll(trackList);
    }

    public void addAll(int index, Collection<Track> trackList) {
        tracks.addAll(index, trackList);
    }

    public void remove(Track track) {
        tracks.remove(track);
    }

    /**
     * Return a composite score for the entire group.  For now use the maximum track
     * score.   Note that scores for tracks not appropriate to the score type will
     * return -Float.MAX, so they are effectively ignored.
     *
     * @param chr
     * @param start
     * @param end
     * @param zoom
     * @param type
     * @param frameName
     * @return
     */
    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, String frameName) {
        float score = -Float.MAX_VALUE;
        for (Track track : tracks) {
            if (track.isVisible()) {
                score = Math.max(score, track.getRegionScore(chr, start, end, zoom, type, frameName));

            }
        }
        return score;
    }


    public String getName() {
        return name;
    }

    public boolean isVisible() {
       return tracks.stream().anyMatch(t -> t != null && t.isVisible());
    }

    public int getHeight() {

        int h = 0;
        for (Track track : tracks) {
            if (track != null && track.isVisible()) {
                h += track.getContentHeight();
            }
        }
        return h;
    }




}