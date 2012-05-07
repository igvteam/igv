/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
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
