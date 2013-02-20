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

/*
 * Created by JFormDesigner on Thu Jun 14 08:42:46 EDT 2012
 */

package org.broad.igv.ui.panel;

import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;

import java.awt.*;
import java.util.Collection;

/**
 * Dialog used for selecting a single FeatureTrack
 * @author Jacob Silterra
 **/
public class FeatureTrackSelectionDialog extends TrackSelectionDialog {

    public FeatureTrackSelectionDialog(Frame owner) {
        super(owner, SelectionMode.SINGLE, IGV.getInstance().getFeatureTracks());
    }

    /**
     * Get the selected track, or null if none
     * @return
     */
    public FeatureTrack getSelectedTrack() {
        Collection<Track> tracks = getSelectedTracks();
        if(tracks == null || tracks.size() == 0) return null;
        return (FeatureTrack) tracks.iterator().next();
    }
}
