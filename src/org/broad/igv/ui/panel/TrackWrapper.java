/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.ui.panel;

import org.broad.igv.track.FeatureTrack;

import java.util.ArrayList;
import java.util.List;

public class TrackWrapper {
    private FeatureTrack track;

    public TrackWrapper(FeatureTrack track) {
        this.track = track;
    }

    public String toString() {
        return track.getName();
    }

    public FeatureTrack getTrack() {
        return this.track;
    }

    public static List<TrackWrapper> wrapTracks(Iterable<FeatureTrack> tracks) {
        ArrayList<TrackWrapper> wrappers = new ArrayList<TrackWrapper>();
        for (FeatureTrack t : tracks) {
            TrackWrapper trackWrapper = new TrackWrapper(t);
            wrappers.add(trackWrapper);
        }
        return wrappers;
    }

}