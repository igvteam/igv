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

import org.apache.log4j.Logger;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackGroup;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.panel.DataPanelContainer;
import org.broad.igv.ui.panel.TrackPanel;
import org.broad.igv.ui.panel.TrackPanelScrollPane;

import java.awt.event.ActionEvent;
import java.util.Collection;
import java.util.List;

/**
 * @author jrobinso
 */
public class FitDataToWindowMenuAction extends MenuAction {

    static Logger log = Logger.getLogger(FitDataToWindowMenuAction.class);
    IGV mainFrame;

    public FitDataToWindowMenuAction(String label, int mnemonic, IGV mainFrame) {
        super(label, null, mnemonic);
        this.mainFrame = mainFrame;
    }

    @Override
    /**
     * The action method. A swing worker is used, so "invoke later" and explicit
     * threads are not neccessary.
     *
     */
    public void actionPerformed(ActionEvent e) {

        for (TrackPanel tp : IGV.getInstance().getTrackPanels()) {
            fitTracksToPanel(tp.getScrollPane().getDataPanel());
        }
        mainFrame.doRefresh();

    }

    /**
     * Adjust the height of  tracks so that all tracks fit in the available
     * height of the panel.  This is not possible in all cases as the
     * minimum height for tracks is respected.
     *
     * @param dataPanel
     * @return
     */
    private boolean fitTracksToPanel(DataPanelContainer dataPanel) {

        boolean success = true;

        int availableHeight = dataPanel.getVisibleHeight();
        int visibleTrackCount = 0;

        // Process data tracks first
        Collection<TrackGroup> groups = dataPanel.getTrackGroups();


        // Count visible tracks.
        for (TrackGroup group : groups) {
            List<Track> tracks = group.getTracks();
            for (Track track : tracks) {
                if (track.isVisible()) {
                    ++visibleTrackCount;
                }
            }
        }


        // Auto resize the height of the visible tracks
        if (visibleTrackCount > 0) {
            int groupGapHeight = (groups.size() + 1) * UIConstants.groupGap;
            double adjustedAvailableHeight = Math.max(1, availableHeight - groupGapHeight);

            double delta = adjustedAvailableHeight / visibleTrackCount;

            // Minimum track height is 1
            if (delta < 1) {
                delta = 1;
            }

            int iTotal = 0;
            double target = 0;
            for (TrackGroup group : groups) {
                List<Track> tracks = group.getTracks();
                for (Track track : tracks) {
                    target += delta;
                    int height = (int) (target - iTotal);
                    track.setHeight(height);
                    iTotal += height;
                }
            }

        }

        return success;
    }

}
