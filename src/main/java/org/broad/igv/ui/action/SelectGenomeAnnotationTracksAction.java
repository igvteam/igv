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
package org.broad.igv.ui.action;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.genome.load.TrackConfig;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.track.Track;
import org.broad.igv.ucsc.Hub;
import org.broad.igv.ucsc.TrackConfigGroup;
import org.broad.igv.ucsc.HubTrackSelectionDialog;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ResourceLocator;

import java.awt.event.ActionEvent;
import java.util.*;

/**
 * @author jrobinso
 */
public class SelectGenomeAnnotationTracksAction extends MenuAction {

    static Logger log = LogManager.getLogger(SelectGenomeAnnotationTracksAction.class);
    IGV mainFrame;

    // Keep track of authorization failures so user isn't constantly harranged
    static HashSet<String> failedURLs = new HashSet();


    public SelectGenomeAnnotationTracksAction(String label, IGV mainFrame) {
        super(label, null);
        this.mainFrame = mainFrame;
    }

    @Override
    public void actionPerformed(ActionEvent evt) {

        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        Hub hub = genome.getHub();
        if (hub == null) {
            // This should not happen
            MessageUtils.showMessage("No annotation tracks available for current genome.");
        }

        final List<Track> loadedTracks = IGV.getInstance().getAllTracks().stream().filter(t -> t.getResourceLocator() != null).toList();
        Set<String> loadedTrackPaths = new HashSet<>(loadedTracks.stream().map(t -> t.getResourceLocator().getPath()).toList());
        List<TrackConfigGroup> groups = hub.getGroupedTrackConfigurations();
        for (TrackConfigGroup g : groups) {
            for (TrackConfig config : g.tracks) {
                config.visible = loadedTrackPaths.contains(config.url);
            }
        }

        HubTrackSelectionDialog dlg = new HubTrackSelectionDialog(groups, IGV.getInstance().getMainFrame());
        dlg.setVisible(true);

        // The dialog action will modify the visible state for each track config
        Set<String> trackPathsToRemove = new HashSet<>();
        List<TrackConfig> tracksToLoad = new ArrayList<>();
        List<TrackConfig> selected = new ArrayList<>();
        for (TrackConfigGroup g : groups) {
            for (TrackConfig config : g.tracks) {
                if (config.visible) {
                    selected.add(config);
                    if (!loadedTrackPaths.contains(config.url)) {
                        tracksToLoad.add(config);
                    }
                } else {
                    trackPathsToRemove.add(config.url);
                }
            }
        }

        List<Track> tracksToRemove = loadedTracks.stream().filter(t -> trackPathsToRemove.contains(t.getResourceLocator().getPath())).toList();
        IGV.getInstance().deleteTracks(tracksToRemove);

        List<ResourceLocator> locators = tracksToLoad.stream().map(t -> ResourceLocator.fromTrackConfig(t)).toList();
        IGV.getInstance().loadTracks(locators);

        // Update genome
        genome.setAnnotationResources(locators);

        // Update preferences
        String key = "hub:" + hub.getUrl();
        PreferencesManager.getPreferences().put(key, String.join(",", selected.stream().map(c -> c.name).toList()));

    }


}
