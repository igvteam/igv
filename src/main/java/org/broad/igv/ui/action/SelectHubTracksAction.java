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
import org.broad.igv.ucsc.hub.TrackConfigContainer;
import org.broad.igv.ucsc.hub.TrackHubSelectionDialog;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.util.*;

/**
 * Select tracks from a track hub.
 *
 * @author jrobinso
 */
public class SelectHubTracksAction extends MenuAction {

    static Logger log = LogManager.getLogger(SelectHubTracksAction.class);

    private Hub hub;
    IGV mainFrame;


    public SelectHubTracksAction(String label, IGV mainFrame, Hub hub) {
        super(label, null);
        this.mainFrame = mainFrame;
        this.hub = hub;
    }

    @Override
    public void actionPerformed(ActionEvent evt) {

        final WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();

        SwingWorker worker = new SwingWorker() {

            @Override
            protected Object doInBackground() throws Exception {
                try {
                    selectTracks(hub);

                } catch (Exception e) {
                    log.error("Error loading track configurations", e);
                    MessageUtils.showMessage(e.getMessage());
                } finally {
                    if (token != null) {
                        WaitCursorManager.removeWaitCursor(token);
                    }
                }
                return null;
            }

            @Override
            protected void done() {
                WaitCursorManager.removeWaitCursor(token);
            }
        };

        worker.execute();

    }

    public static void selectTracks(Hub hub) {

        final List<Track> loadedTracks = IGV.getInstance().getAllTracks().stream().filter(t -> t.getResourceLocator() != null).toList();
        Set<String> loadedTrackPaths = new HashSet<>(loadedTracks.stream().map(t -> t.getResourceLocator().getPath()).toList());

        TrackHubSelectionDialog dlg = TrackHubSelectionDialog.getTrackHubSelectionDialog(hub, loadedTrackPaths, false);

        dlg.setVisible(true);

        if (!dlg.isCanceled()) {

            List<TrackConfigContainer> groups = hub.getGroupedTrackConfigurations();

            // The dialog action will modify the visible state for each track config
            List<TrackConfig> tracksToLoad = new ArrayList<>();
            List<TrackConfig> selected = new ArrayList<>();
            for (TrackConfigContainer g : groups) {
                g.map(config -> {
                    if (config.getVisible()) {
                        selected.add(config);
                        if (!loadedTrackPaths.contains(config.getUrl())) {
                            tracksToLoad.add(config);
                        }
                    }
                    return null;
                });
            }

            List<ResourceLocator> locators = tracksToLoad.stream().map(t -> ResourceLocator.fromTrackConfig(t)).toList();
            IGV.getInstance().loadTracks(locators);
        }
    }

}
