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

package org.broad.igv.ui;

import org.broad.igv.track.Track;
import org.broad.igv.util.Filter;
import org.broad.igv.util.FilterElement;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * @author eflakes
 */
public class TrackFilter extends Filter {

    private boolean showAll = false;

    List<Track> currentTrackList;

    public TrackFilter(String name, List<TrackFilterElement> elements) {
        super(name, elements);
    }

    public void clearTracks() {
        currentTrackList = null;
    }

    public void addTracks(List<Track> tracks) {

        if (tracks == null || tracks.isEmpty())
            return;

        if (currentTrackList == null) {
            currentTrackList = new ArrayList<Track>(tracks);
        }
    }

    /**
     * Evaluate the TrackFilterElement set.
     *
     * @return
     */
    @Override
    public void evaluate() {

        if (currentTrackList == null || currentTrackList.isEmpty())
            return;

        boolean filterEnabled = isEnabled();

        for (Track track : currentTrackList) {

            // Must start as null which means no previous results
            Boolean result = null;

            // If filter is not enabled just show all tracks
            if (!filterEnabled) {
                track.setVisible(!filterEnabled);
                continue;
            }

            if (track.isFilterable()) {

                // Evaluate tracks
                Iterator iterator = getFilterElements();
                while (iterator.hasNext()) {

                    FilterElement element = (FilterElement) iterator.next();
                    result = element.evaluate(track, result);
                }
                track.setVisible(result);
            }
        }
    }


}
