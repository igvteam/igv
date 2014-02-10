package org.broad.igv.cursor;

import java.util.*;

/**
 * The "model" object for a Cursor instance
 *
 * @author jrobinso
 *         Date: 1/14/14
 *         Time: 12:43 PM
 */
public class CursorModel {

    public static int frameBPWidth = 1000;

    private List<CursorTrack> tracks;
    private List<CursorRegion> regions;
    private List<CursorRegion> filteredRegions;
    private RegionFilter filter;
    private double framePixelWidth = 24;
    int frameMargin = 6;
    private double origin = 0;
    private int framePixelHeight = 50;
    CursorTrack sortedTrack;

    public List<CursorTrack> getTracks() {
        return tracks;
    }

    public void setTracks(List<CursorTrack> tracks) {
        this.tracks = tracks;
    }

    public void setRegions(List<CursorRegion> frames) {
        this.regions = frames;
        updateFilteredRegions();
    }

    private void updateFilteredRegions() {
        if(filter == null || regions == null) filteredRegions = null;
        else {
            filteredRegions = new ArrayList<CursorRegion>();
            for(CursorRegion r : regions) {
                if(filter.pass(r)) filteredRegions.add(r);
            }
        }
    }

    public List<CursorRegion> getFilteredRegions() {
        return filteredRegions == null ? regions : filteredRegions;
    }

    public RegionFilter getFilter() {
        return filter;
    }

    public void setFilter(RegionFilter filter) {
        this.filter = filter;
        updateFilteredRegions();
    }

    public double getFramePixelWidth() {
        return framePixelWidth;
    }

    public void setFramePixelWidth(double framePixelWidth) {
        this.framePixelWidth = framePixelWidth;
        this.frameMargin = (int) Math.min(8, framePixelWidth / 4);
    }

    public int getFrameBPWidth() {
        return frameBPWidth;
    }

    public void setFrameBPWidth(int frameBPWidth) {
        this.frameBPWidth = frameBPWidth;
    }

    public double getOrigin() {
        return origin;
    }

    public void setOrigin(double origin) {
        this.origin = origin;
    }

    public int getTrackPixelHeight() {
        return framePixelHeight;
    }

    public void setFramePixelHeight(int framePixelHeight) {
        this.framePixelHeight = framePixelHeight;
    }

    public void addTrack(CursorTrack t) {
        if(tracks == null) tracks = new ArrayList<CursorTrack>();
        tracks.add(t);
    }

    // Sort frames based on signal from track t
    public void sortFrames(final CursorTrack t, final int sortDirection) {

        sortedTrack = t;
        // First, randomize the frames to prevent memory from previous sorts.  There are many ties (e.g. zeroes)
        // so a stable sort carries a lot of memory, which can be confusing and imply correlations where none exist.
        Collections.shuffle(regions);
        Collections.sort(regions, new Comparator<CursorRegion>() {

            @Override
            public int compare(CursorRegion cursorRegion1, CursorRegion cursorRegion2) {

                int l1 = t.getLongestFeatureLength(cursorRegion1.getChr());
                int l2 = t.getLongestFeatureLength(cursorRegion2.getChr());
                double s1 = cursorRegion1.getScore(t.getFeatures(cursorRegion1.getChr()), l1, frameBPWidth);
                double s2 = cursorRegion2.getScore(t.getFeatures(cursorRegion2.getChr()), l2, frameBPWidth);
                return sortDirection * (s1 == s2 ? 0 : (s1 > s2 ? -1 : 1));

            }
        });

    }

    public void shiftOriginPixels(int delta) {

        origin = Math.max(0, origin + ((double) delta) / framePixelWidth);

    }

    public int getFrameMargin() {
        return frameMargin;
    }

    public CursorTrack getSortedTrack() {
        return sortedTrack;
    }

}
