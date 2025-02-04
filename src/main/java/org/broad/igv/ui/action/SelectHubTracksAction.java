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

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.genome.load.TrackConfig;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.track.Track;
import org.broad.igv.ucsc.hub.Hub;
import org.broad.igv.ucsc.hub.TrackConfigGroup;
import org.broad.igv.ucsc.hub.TrackHubSelectionDialog;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ResourceLocator;

import java.awt.event.ActionEvent;
import java.util.*;

/**
 * Select tracks from a track hub.  This action is used in 2 modes, (1) load tracks from a specific hub, and (2) select
 * default annotation tracks for the currently loaded genome.  Tracks are loaded in both modes, but in the
 * second selected tracks are also added to the genome definition.
 *
 * @author jrobinso
 */
public class SelectHubTracksAction extends MenuAction {

    static Logger log = LogManager.getLogger(SelectHubTracksAction.class);
    private Hub hub;
    IGV mainFrame;

    // Keep track of authorization failures so user isn't constantly harranged
    static HashSet<String> failedURLs = new HashSet();

    boolean updateGenome;

    public SelectHubTracksAction(String label, IGV mainFrame, Hub hub) {
        super(label, null);
        this.mainFrame = mainFrame;
        this.updateGenome = hub == null;
        this.hub = hub;
    }

    @Override
    public void actionPerformed(ActionEvent evt) {

        try {
            Genome genome = GenomeManager.getInstance().getCurrentGenome();
            if (hub == null) {
                hub = genome.getGenomeHub();
                if (hub == null) {
                    // This should not happen
                    MessageUtils.showMessage("No tracks available for current genome.");
                }
            }

            final List<Track> loadedTracks = IGV.getInstance().getAllTracks().stream().filter(t -> t.getResourceLocator() != null).toList();
            Set<String> loadedTrackPaths = new HashSet<>(loadedTracks.stream().map(t -> t.getResourceLocator().getPath()).toList());
            List<TrackConfigGroup> groups = hub.getGroupedTrackConfigurations();
            for (TrackConfigGroup g : groups) {
                for (TrackConfig config : g.tracks) {
                    config.setVisible(loadedTrackPaths.contains(config.getUrl()));
                }
            }

            TrackHubSelectionDialog dlg = new TrackHubSelectionDialog(hub, groups, IGV.getInstance().getMainFrame());
            dlg.setVisible(true);

            // The dialog action will modify the visible state for each track config
            Set<String> trackPathsToRemove = new HashSet<>();
            List<TrackConfig> tracksToLoad = new ArrayList<>();
            List<TrackConfig> selected = new ArrayList<>();
            for (TrackConfigGroup g : groups) {
                for (TrackConfig config : g.tracks) {
                    if (config.getVisible()) {
                        selected.add(config);
                        if (!loadedTrackPaths.contains(config.getUrl())) {
                            tracksToLoad.add(config);
                        }
                    } else {
                        trackPathsToRemove.add(config.getUrl());
                    }
                }
            }

            List<Track> tracksToRemove = loadedTracks.stream().filter(t -> trackPathsToRemove.contains(t.getResourceLocator().getPath())).toList();
            IGV.getInstance().deleteTracks(tracksToRemove);

            List<ResourceLocator> locators = tracksToLoad.stream().map(t -> ResourceLocator.fromTrackConfig(t)).toList();
            IGV.getInstance().loadTracks(locators);

            // Update genome
            if (updateGenome) {
                genome.setAnnotationResources(locators);
                // Update preferences
                String key = "hub:" + hub.getUrl();
                PreferencesManager.getPreferences().put(key, String.join(",", selected.stream().map(c -> c.getName()).toList()));
            }
        } catch (Exception e) {
            log.error("Error loading track configurations", e);
            MessageUtils.showMessage(e.getMessage());
        }

    }


}
