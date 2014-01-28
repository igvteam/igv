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

    private List<CursorTrack> tracks;
    private List<CursorRegion> cursorRegions;
    private int framePixelWidth = 24;
    int frameMargin = 6;
    private int frameBPWidth = 1000;
    private double origin = 0;
    private int framePixelHeight = 50;
    CursorTrack sortedTrack;

    public List<CursorTrack> getTracks() {
        return tracks;
    }

    public void setTracks(List<CursorTrack> tracks) {
        this.tracks = tracks;
    }

    public List<CursorRegion> getFrames() {
        return cursorRegions;
    }

    public void setFrames(List<CursorRegion> frames) {
        this.cursorRegions = frames;
    }

    public int getFramePixelWidth() {
        return framePixelWidth;
    }

    public void setFramePixelWidth(int framePixelWidth) {
        this.framePixelWidth = framePixelWidth;
        this.frameMargin = Math.min(8, Math.max(1, framePixelWidth / 4));
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
        Collections.shuffle(cursorRegions);

        Collections.sort(cursorRegions, new Comparator<CursorRegion>() {

            @Override
            public int compare(CursorRegion cursorRegion1, CursorRegion cursorRegion2) {
                double s1 = cursorRegion1.getScore(t.getFeatures(cursorRegion1.getChr()));
                double s2 = cursorRegion2.getScore(t.getFeatures(cursorRegion2.getChr()));
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
