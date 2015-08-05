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

package org.broad.igv.cursor;

import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.SignalFeature;

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
    CursorMainWindow mainWindow;

    public CursorModel(CursorMainWindow cursorMainWindow) {
        this.mainWindow = cursorMainWindow;
    }

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
        if (filter == null || regions == null) filteredRegions = null;
        else {
            filteredRegions = new ArrayList<CursorRegion>();
            for (CursorRegion r : regions) {
                if (filter.pass(r)) filteredRegions.add(r);
            }
        }
        mainWindow.updateRegionsLabel();
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
        mainWindow.updateRegionsLabel();
    }

    public int getFrameBPWidth() {
        return frameBPWidth;
    }

    public void setFrameBPWidth(int frameBPWidth) {
        this.frameBPWidth = frameBPWidth;
        mainWindow.updateRegionsLabel();
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
        if (tracks == null) tracks = new ArrayList<CursorTrack>();
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

                double s1 = cursorRegion1.getScore(t, frameBPWidth);
                double s2 = cursorRegion2.getScore(t, frameBPWidth);
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
