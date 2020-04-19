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

    //static Logger log = Logger.getLogger(GroupTracksMenuAction.class);
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
                IGV.getInstance().removeTracks(Arrays.asList(mergedTracks));
            }
        }
        IGV.getInstance().repaint();
    }

}
