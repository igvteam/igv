/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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
package org.broad.igv.ui.panel;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackClickEvent;
import org.broad.igv.track.TrackMenuUtils;
import org.broad.igv.ui.IGVMainFrame;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * @author eflakes
 */
abstract public class TrackPanelComponent extends JPanel {

    private static Logger log = Logger.getLogger(TrackPanelComponent.class);
    List<MouseableRegion> trackRegions;

    private TrackPanel trackPanel;


    public TrackPanelComponent(TrackPanel trackPanel) {
        this.trackPanel = trackPanel;
        setFocusable(true);
        trackRegions = new ArrayList();
    }

    public TrackPanel getTrackSetView() {
        if (trackPanel == null) {
            trackPanel = (TrackPanel) getParent();
        }
        return trackPanel;
    }

    public String getTrackSetID() {
        return getTrackSetView().getName();
    }


    protected void addMousableRegion(MouseableRegion region) {
        trackRegions.add(region);
    }

    protected void removeMousableRegions() {
        trackRegions.clear();
    }

    protected List<MouseableRegion> getTrackRegions() {
        return trackRegions;
    }

    public boolean scrollTo(String trackName) {
        Track t = findNextTrackMatching(trackName);
        if (t != null) {
            IGVMainFrame.getInstance().getTrackManager().clearSelections();
            t.setSelected(true);
            if (trackPanel.getScrollPane().getVerticalScrollBar().isShowing()) {
                trackPanel.getScrollPane().getVerticalScrollBar().setValue(t.getTop());
            }
            return true;
        }

        return false;
    }

    int searchIdx = 0;

    private synchronized Track findNextTrackMatching(String trackName) {
        List<Track> tracks = getAllTracks();
        searchIdx = Math.min(searchIdx, tracks.size());
        for (int i = searchIdx; i < tracks.size(); i++) {
            Track t = tracks.get(i);
            if (t.getName().toUpperCase().contains(trackName.toUpperCase())) {
                searchIdx = i + 1;
                return t;
            }
        }
        for (int i = 0; i < searchIdx; i++) {
            Track t = tracks.get(i);
            if (t.getName().toUpperCase().contains(trackName.toUpperCase())) {
                searchIdx = i + 1;
                return t;
            }
        }
        return null;
    }

    /**
     * Method description
     *
     * @param x
     * @param y
     * @return
     */
    abstract public String getPopupMenuTitle(int x, int y);

    protected Collection<Track> getSelectedTracks() {
        return IGVMainFrame.getInstance().getTrackManager().getSelectedTracks();
    }

    public List<Track> getAllTracks() {
        TrackPanel dataTrackView = (TrackPanel) getParent();
        return dataTrackView.getTracks();
    }


    protected void openPopupMenu(TrackClickEvent te) {
        MouseEvent e = te.getMouseEvent();


        // If there is a single selected track give it a chance to handle the click
        Collection<Track> selectedTracs = getSelectedTracks();

        if (selectedTracs.size() == 0) {
            return;
        }

        if (selectedTracs.size() == 1) {
            if (selectedTracs.iterator().next().handleClick(te)) {
                return;
            }
        }
        String title = getPopupMenuTitle(e.getX(), e.getY());

        JPopupMenu menu = TrackMenuUtils.getPopupMenu(selectedTracs, title, te);

        menu.show(e.getComponent(), e.getX(), e.getY());
    }

    protected void toggleTrackSelections(MouseEvent e) {
        for (MouseableRegion mouseRegion : trackRegions) {
            if (mouseRegion.containsPoint(e.getX(), e.getY())) {
                IGVMainFrame.getInstance().getTrackManager().toggleTrackSelections(mouseRegion.getTracks());
                return;
            }
        }
    }

    protected void clearTrackSelections() {
        IGVMainFrame.getInstance().getTrackManager().clearSelections();
    }

    protected void selectTracks(MouseEvent e) {

        if (log.isDebugEnabled()) {
            log.debug("Enter selectTracks");
        }



        for (MouseableRegion mouseRegion : trackRegions) {
            if (mouseRegion.containsPoint(e.getX(), e.getY())) {
                IGVMainFrame.getInstance().getTrackManager().setTrackSelections(mouseRegion.getTracks());
                return;
            }
        }
    }

    protected boolean isTrackSelected(MouseEvent e) {
        for (MouseableRegion mouseRegion : trackRegions) {
            if (mouseRegion.containsPoint(e.getX(), e.getY())) {
                for (Track t : mouseRegion.getTracks()) {
                    if (t.isSelected()) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

}
