package org.igv.ui.action;

import org.igv.feature.genome.load.TrackConfig;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.ucsc.hub.*;
import org.igv.ui.IGV;
import org.igv.ui.IGVMenuBar;
import org.igv.ui.WaitCursorManager;
import org.igv.ui.util.MessageUtils;
import org.igv.ui.util.UIUtilities;
import org.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Select tracks from a track hub.
 *
 * @author jrobinso
 */
public class SelectHubTracksAction extends MenuAction {

    static Logger log = LogManager.getLogger(SelectHubTracksAction.class);

    private HubDescriptor hubDescriptor;
    String genomeId;

    public SelectHubTracksAction(String label, HubDescriptor hubDescriptor, String id) {
        super(label, null);
        this.hubDescriptor = hubDescriptor;
        this.genomeId = id;
    }

    @Override
    public void actionPerformed(ActionEvent evt) {

        final WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();

        SwingWorker worker = new SwingWorker() {

            @Override
            protected Object doInBackground() throws Exception {
                try {
                    Hub hub = HubParser.loadHub(hubDescriptor.getUrl());
                    if(hub.getSupportedTrackCount(genomeId) > 0) {
                        selectAndLoadTracks(hub, genomeId);
                    } else {
                        MessageUtils.showMessage(hubDescriptor.getShortLabel() +
                                " does not have any IGV supported tracks for " + genomeId);
                    }

                } catch (Exception e) {
                    log.error("Error selecting or loading tracks for hub " + hubDescriptor.getShortLabel(), e);
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

    public static void selectAndLoadTracks(Hub hub, String id) {

        Set<String> loadedTrackPaths = IGV.getInstance().getAllTracks().stream()
                .filter(t -> t.getResourceLocator() != null)
                .map(t -> t.getResourceLocator().getPath())
                .collect(Collectors.toSet());

        // Called from background threads; the dialog is built and shown on the event thread.
        final TrackSelectionDialog[] holder = new TrackSelectionDialog[1];
        UIUtilities.invokeAndWaitOnEventThread(() -> {
            holder[0] = TrackSelectionDialog.getTrackHubSelectionDialog(hub, id, loadedTrackPaths, false, null);
            holder[0].setVisible(true);
        });
        TrackSelectionDialog dlg = holder[0];

        if (!dlg.isCanceled()) {

            List<TrackConfigContainer> groups = hub.getGroupedTrackConfigurations(id);

            // The dialog action will modify the visible state for each track config.  Checked tracks not yet loaded
            // are loaded, loaded tracks that were unchecked are unloaded.
            List<TrackConfig> tracksToLoad = new ArrayList<>();
            Set<String> selectedPaths = new HashSet<>();
            Set<String> trackPathsToRemove = new HashSet<>();
            for (TrackConfigContainer g : groups) {
                g.map(config -> {
                    if (config.visible) {
                        selectedPaths.add(config.url);
                        if (!loadedTrackPaths.contains(config.url)) {
                            tracksToLoad.add(config);
                        }
                    } else if (loadedTrackPaths.contains(config.url)) {
                        trackPathsToRemove.add(config.url);
                    }
                    return null;
                });
            }
            // Entries can share a URL; keep it if any entry using it is still checked
            trackPathsToRemove.removeAll(selectedPaths);

            IGV.getInstance().deleteTracksByPath(trackPathsToRemove);

            List<ResourceLocator> locators = tracksToLoad.stream().map(t -> ResourceLocator.fromTrackConfig(t)).toList();
            IGV.getInstance().loadTracks(locators);
        }
    }

}
