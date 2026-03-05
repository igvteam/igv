package org.igv.ui.panel;


import org.igv.Globals;
import org.igv.track.Track;
import org.igv.track.TrackClickEvent;
import org.igv.track.TrackMenuUtils;
import org.igv.ui.IGV;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;


/**
 * @author eflakes
 */
abstract public class TrackPanelComponent extends JPanel {

    protected final boolean darkMode;
    List<MouseableRegion> mouseRegions;
    private TrackPanel trackPanel;

    public TrackPanelComponent(TrackPanel trackPanel) {
        this.trackPanel = trackPanel;
        setFocusable(true);
        mouseRegions = new ArrayList<>();
        this.darkMode = Globals.isDarkMode();
    }

    public TrackPanel getTrackPanel() {
        return trackPanel;
    }

    public Track getTrack() {
        return getTrackPanel().getTrack();
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

    public List<Track> getAllTracks() {
        TrackPanel dataTrackView = (TrackPanel) getParent();
        return dataTrackView.getTracks();
    }

    protected void openPopupMenu(TrackClickEvent te) {

        MouseEvent e = te.getMouseEvent();
        IGVPopupMenu menu = null;

        Track track = getTrack();
        if (track == null) {
            // If this panel is empty (no tracks), and is removable, present option to remove it
            if (getAllTracks().size() == 0) {
                menu = new IGVPopupMenu();
                JMenuItem item = new JMenuItem("Remove panel");
                item.addActionListener(e12 -> {
                    IGV.getInstance().removeTrackPanel(getTrackPanel());
                });
                menu.add(item);
            }
        } else {
            Collection<Track> tracks = Collections.singletonList(track);

            // Give the track an opportunity to provide the popup menu
            menu = track.getPopupMenu(te);

            // If still no menu, create a generic one with common items
            if (menu == null) {
                String title = track.getName();
                menu = TrackMenuUtils.getPopupMenu(tracks, title, te);
            }

            // Add additional items, if any
            if (menu.includeStandardItems()) {

                TrackMenuUtils.addPluginItems(menu, tracks, te);

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
                JMenuItem exportFeats = TrackMenuUtils.getExportFeatures(tracks, frame);
                if (exportFeats != null) menu.add(exportFeats);

                menu.addSeparator();
                menu.add(TrackMenuUtils.getRemoveMenuItem(tracks));
            }
        }

        if (menu != null) {
            menu.show(e.getComponent(), e.getX(), e.getY());
        }

    }


    public void saveImage(String extension) {
        IGV.getInstance().saveImage(getTrackPanel().getScrollPane(), "igv_panel", extension);
    }


}
