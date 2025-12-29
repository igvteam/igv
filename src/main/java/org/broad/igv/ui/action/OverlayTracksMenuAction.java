/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.action;

import org.broad.igv.track.DataTrack;
import org.broad.igv.track.MergedTracks;
import org.broad.igv.track.Track;
import org.broad.igv.ui.AttributeSelectionDialog;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.TrackPanel;
import org.broad.igv.ui.util.UIUtilities;

import java.awt.event.ActionEvent;
import java.util.*;

/**
 * @author jrobinso
 */
public class OverlayTracksMenuAction extends MenuAction {

    //static Logger log = LogManager.getLogger(GroupTracksMenuAction.class);
    IGV igv;

    public OverlayTracksMenuAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
    }

    @Override
    public void actionPerformed(ActionEvent e) {

        UIUtilities.invokeOnEventThread(() -> {

            final AttributeSelectionDialog dlg = new AttributeSelectionDialog(igv.getMainFrame(), "Overlay");
            dlg.setVisible(true);

            if (!dlg.isCanceled()) {
                String selectedAttribute = dlg.getSelected();
                if (selectedAttribute == null) {
                    unmerge(IGV.getInstance().getAllTracks());
                } else {
                    List<DataTrack> tracks = IGV.getInstance().getDataTracks();
                    Map<String, List<DataTrack>> groups = new HashMap<>();
                    for (DataTrack t : tracks) {
                        String v = t.getAttributeValue(selectedAttribute);
                        if (v != null) {
                            List<DataTrack> tlist = groups.get(v);
                            if (tlist == null) {
                                tlist = new ArrayList<>();
                                groups.put(v, tlist);
                            }
                            tlist.add(t);
                        }
                    }

                    for (Map.Entry<String, List<DataTrack>> entry : groups.entrySet()) {
                        String name = entry.getKey();
                        merge(entry.getValue(), name);
                    }
                    igv.repaint();
                }

            }
        });
    }

    public static void merge(List<DataTrack> dataTrackList, String name) {
        MergedTracks mergedTracks = new MergedTracks(UUID.randomUUID().toString(), name, dataTrackList);
        Track firstTrack = dataTrackList.iterator().next();
        TrackPanel panel = TrackPanel.getParentPanel(firstTrack);
        panel.addTrack(mergedTracks);
        panel.moveSelectedTracksTo(Arrays.asList(mergedTracks), firstTrack, false);
        panel.removeTracks(dataTrackList);
    }

    public static void unmerge(Collection<Track> tracks) {
        for (Track t : tracks) {
            if (t instanceof MergedTracks) {
                TrackPanel panel = TrackPanel.getParentPanel(t);
                MergedTracks mergedTracks = (MergedTracks) t;
                mergedTracks.setTrackAlphas(1.0);
                panel.addTracks(mergedTracks.getMemberTracks());
                panel.moveSelectedTracksTo(mergedTracks.getMemberTracks(), mergedTracks, true);
                IGV.getInstance().deleteTracks(Arrays.asList(mergedTracks));
            }
        }
        IGV.getInstance().repaint();
    }

}
