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
package org.broad.igv.ui.panel;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.logging.*;
import org.broad.igv.Globals;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackClickEvent;
import org.broad.igv.track.TrackMenuUtils;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;


/**
 * @author eflakes
 */
abstract public class TrackPanelComponent extends JPanel {

    private static final String DELETE_TRACKS_KEY = "deleteTracks";

    protected final boolean darkMode;
    List<MouseableRegion> mouseRegions;
    private TrackPanel trackPanel;

    public TrackPanelComponent(TrackPanel trackPanel) {
        this.trackPanel = trackPanel;
        setFocusable(true);
        mouseRegions = new ArrayList<>();

        initKeyDispatcher();
        this.darkMode = Globals.isDarkMode();
    }

    private void initKeyDispatcher() {
        final Action delTracksAction = new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent e) {
                TrackMenuUtils.removeTracksAction(IGV.getInstance().getSelectedTracks());
            }
        };

        final KeyStroke delKey = KeyStroke.getKeyStroke(KeyEvent.VK_DELETE, 0, false);
        final KeyStroke backspaceKey = KeyStroke.getKeyStroke(KeyEvent.VK_BACK_SPACE, 0, false);

        getInputMap().put(delKey, DELETE_TRACKS_KEY);
        getInputMap().put(backspaceKey, DELETE_TRACKS_KEY);
        getActionMap().put(DELETE_TRACKS_KEY, delTracksAction);

    }

    public TrackPanel getTrackPanel() {
        if (trackPanel == null) {
            trackPanel = (TrackPanel) getParent();
        }
        return trackPanel;
    }

    public String getTrackSetID() {
        return getTrackPanel().getName();
    }


    protected void addMousableRegion(MouseableRegion region) {
        mouseRegions.add(region);
    }

    protected void removeMousableRegions() {
        mouseRegions.clear();
    }

    protected List<MouseableRegion> getMouseRegions() {
        return mouseRegions;
    }

    public boolean scrollTo(String trackName) {
        Track t = findNextTrackMatching(trackName);
        if (t != null) {
            IGV.getInstance().clearSelections();
            t.setSelected(true);
            scrollToPosition(t.getY());
            return true;
        }
        return false;
    }

    public void scrollToPosition(int y) {
        if (trackPanel.getScrollPane().getVerticalScrollBar().isShowing()) {
            trackPanel.getScrollPane().getVerticalScrollBar().setValue(y);
        }
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

    public String getPopupMenuTitle(int x, int y) {

        Collection<Track> tracks = getSelectedTracks();

        String popupTitle;
        if (tracks.size() == 1) {
            popupTitle = tracks.iterator().next().getName();
        } else {
            popupTitle = "Total Tracks Selected: " + tracks.size();
        }

        return popupTitle;
    }

    protected Collection<Track> getSelectedTracks() {
        return IGV.getInstance().getSelectedTracks();
    }

    public List<Track> getAllTracks() {
        TrackPanel dataTrackView = (TrackPanel) getParent();
        return dataTrackView.getTracks();
    }

    protected void openPopupMenu(TrackClickEvent te) {
        openPopupMenu(te, null);
    }

    protected void openPopupMenu(TrackClickEvent te, List<Component> extraItems) {

        MouseEvent e = te.getMouseEvent();
        IGVPopupMenu menu = null;

        final Collection<Track> selectedTracks = getSelectedTracks();
        if (selectedTracks.isEmpty()) {
            // If this panel is empty (no tracks), and is removable,  present option to remove it
            if (getAllTracks().size() == 0 && IGV.getInstance().panelIsRemovable(getTrackPanel())) {
                menu = new IGVPopupMenu();
                JMenuItem item = new JMenuItem("Remove panel");
                item.addActionListener(e12 -> {
                    IGV.getInstance().removeTrackPanel(getTrackPanel());
                });
                menu.add(item);
            }
        } else {


            // If a single track is selected, give it an opportunity to provide the popup menu
            if (selectedTracks.size() == 1) {
                Track track = selectedTracks.iterator().next();
                menu = track.getPopupMenu(te);
            }

            // If still no menu, create a generic one with common items
            if (menu == null) {
                String title = getPopupMenuTitle(e.getX(), e.getY());
                menu = TrackMenuUtils.getPopupMenu(selectedTracks, title, te);
            }

            // Add additional items, if any
            if (extraItems != null) {
                menu.addSeparator();
                for (Component item : extraItems) {
                    menu.add(item);
                }
            }

            if (menu.includeStandardItems()) {

                TrackMenuUtils.addPluginItems(menu, selectedTracks, te);

                // Add saveImage items
                menu.addSeparator();
                JMenuItem savePng = new JMenuItem("Save PNG image...");
                savePng.addActionListener(e1 -> saveImage("png"));
                menu.add(savePng);
                JMenuItem saveSvg = new JMenuItem("Save SVG image...");
                saveSvg.addActionListener(e1 -> saveImage("svg"));
                menu.add(saveSvg);

                // Add export features
                ReferenceFrame frame = FrameManager.getDefaultFrame();
                JMenuItem exportFeats = TrackMenuUtils.getExportFeatures(selectedTracks, frame);
                if (exportFeats != null) menu.add(exportFeats);

                menu.addSeparator();
                menu.add(TrackMenuUtils.getRemoveMenuItem(selectedTracks));
            }
        }

        if(menu != null) {
            menu.show(e.getComponent(), e.getX(), e.getY());
        }

    }

    protected void toggleTrackSelections(MouseEvent e) {
        for (MouseableRegion mouseRegion : mouseRegions) {
            if (mouseRegion.containsPoint(e.getX(), e.getY())) {
                IGV.getInstance().toggleTrackSelections(mouseRegion.getTracks());
                return;
            }
        }
    }

    protected void clearTrackSelections() {
        IGV.getInstance().clearSelections();
        IGV.getInstance().getMainFrame().repaint();
    }

    protected void selectTracks(MouseEvent e) {
        for (MouseableRegion mouseRegion : mouseRegions) {
            if (mouseRegion.containsPoint(e.getX(), e.getY())) {
                IGV.getInstance().setTrackSelections(mouseRegion.getTracks());
                return;
            }
        }
    }

    protected boolean isTrackSelected(MouseEvent e) {
        for (MouseableRegion mouseRegion : mouseRegions) {
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


    public void saveImage(String extension) {
        IGV.getInstance().saveImage(getTrackPanel().getScrollPane(), "igv_panel", extension);
    }


}
