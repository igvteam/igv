/*
 * Created by JFormDesigner on Thu Jun 14 08:42:46 EDT 2012
 */

package org.broad.igv.ui.panel;

import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;

import java.awt.*;
import java.util.*;

/**
 * Dialog used for selecting a single FeatureTrack
 * @author Jacob Silterra
 **/
public class FeatureTrackSelectionDialog extends TrackSelectionDialog {

    public FeatureTrackSelectionDialog(Frame owner) {
        super(owner, SelectionMode.SINGLE, IGV.getInstance().getFeatureTracks());
    }

    public FeatureTrackSelectionDialog(Frame owner, java.util.List<FeatureTrack> tracks) {
        super(owner, SelectionMode.SINGLE, tracks);
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
