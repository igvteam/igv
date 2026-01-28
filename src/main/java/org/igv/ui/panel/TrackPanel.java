package org.igv.ui.panel;


import org.igv.feature.RegionOfInterest;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.track.RegionScoreType;
import org.igv.track.Track;
import org.igv.track.TrackGroup;
import org.igv.ui.IGV;
import org.igv.ui.UIConstants;

import javax.swing.*;
import java.awt.*;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * @author eflakes
 */
public class TrackPanel extends IGVPanel implements Scrollable {

    private static Logger log = LogManager.getLogger(TrackPanel.class);

    private String name = null;
    //private JPanel grabPanel;
    private TrackNamePanel namePanel;
    private AttributePanel attributePanel;
    private DataPanelContainer dataPanelContainer;
    private List<TrackGroup> trackGroups;
    private Track track;

    transient int lastHeight = 0;

    /**
     * Constructs ...
     *
     * @param name
     */
    public TrackPanel(String name, MainPanel mainPanel) {
        super(mainPanel);
        this.name = name;
        init();
    }


    private void init() {
        setBackground(Color.white);
        namePanel = new TrackNamePanel(this);
        attributePanel = new AttributePanel(this);
        dataPanelContainer = new DataPanelContainer(this);

        //add(grabPanel);
        add(namePanel);
        add(attributePanel);
        add(dataPanelContainer);

    }


    public void createDataPanels() {
        dataPanelContainer.createDataPanels();
    }

    @Override
    public void doLayout() {
        synchronized (getTreeLock()) {
            // Use preferred height (content height) instead of actual height
            // This ensures the panel is sized correctly for scrolling
            int h = getPreferredSize().height;

            int nw = mainPanel.getNamePanelWidth();

            namePanel.setBounds(mainPanel.getNamePanelX(), 0, nw, h);
            attributePanel.setBounds(mainPanel.getAttributePanelX(), 0, mainPanel.getAttributePanelWidth(), h);
            dataPanelContainer.setBounds(mainPanel.getDataPanelX(), 0, mainPanel.getDataPanelWidth(), h);

            dataPanelContainer.doLayout();
        }
    }

    @Override
    public void setBackground(Color color) {
        super.setBackground(color);
        if (namePanel != null) {
            namePanel.setBackground(color);
            attributePanel.setBackground(color);
            dataPanelContainer.setBackground(color);
        }
    }

    public Track getTrack() {
        return track;

    }

    public TrackNamePanel getNamePanel() {
        return namePanel;
    }


    public AttributePanel getAttributePanel() {
        return attributePanel;
    }


    public DataPanelContainer getDataPanelContainer() {
        return dataPanelContainer;
    }


    public String getName() {
        return name;
    }

    /**
     * Method description
     *
     * @return
     */
    public List<TrackGroup> getGroups() {
        return trackGroups;
    }

    /**
     * Method description
     *
     * @return
     */
    public boolean hasTracks() {
        return track != null;
    }

    public int getVisibleTrackCount() {
        return track != null && track.isVisible() ? 1 : 0;
    }

    public boolean isHeightChanged() {
        int height = track.getContentHeight();
        boolean change = height != lastHeight;
        lastHeight = height;
        return change;
    }

    public List<Track> getTracks() {
        return Collections.singletonList(this.track);
    }

    public boolean containsTrack(Track track) {
        return track == this.track;
    }

    public void clearTracks() {
        if (track != null) {
            track.unload();
        }
        track = null;
    }


    public boolean fitTracksToPanel() {
        DataPanelContainer dataPanel = this.getScrollPane().getDataPanel();
        boolean success = true;

        int availableHeight = dataPanel.getVisibleHeight();
        int visibleSampleCount = 0;

        // Process data tracks first
        Collection<TrackGroup> groups = dataPanel.getTrackGroups();


        // Count visible samples.
        for (TrackGroup group : groups) {
            List<Track> tracks = group.getVisibleTracks();
            for (Track track : tracks) {
                if (track.isVisible()) {
                    visibleSampleCount += track.sampleCount();
                }
            }
        }


        // Auto resize the height of the visible tracks
        if (visibleSampleCount > 0) {
            int groupGapHeight = (groups.size() + 1) * UIConstants.groupGap;
            double adjustedAvailableHeight = Math.max(1, availableHeight - groupGapHeight);

            double delta = adjustedAvailableHeight / visibleSampleCount;

            // Minimum track height is 1
            if (delta < 1) {
                delta = 1;
            }

            int iTotal = 0;
            double target = 0;
            for (TrackGroup group : groups) {
                List<Track> tracks = group.getVisibleTracks();
                for (Track track : tracks) {
                    target += delta * track.sampleCount();
                    int height = (int) (target - iTotal);
                    track.setHeight(height);
                    iTotal += height;
                }
            }

        }

        return success;
    }

    /**
     * Add a track to this panel.  If tracks are grouped, search for correct group, or make a new one if not found.
     *
     * @param track
     */
    public void addTrack(Track track) {
        if (this.track != null) {
            throw new RuntimeException("TrackPanel can only hold a single track");
        }
        this.track = track;
    }

    public void addTracks(Collection<? extends Track> tracks) {
        for (Track t : tracks) {
            addTrack(t);
        }
    }

    public void reset() {
        clearTracks();
    }

    public void sortTracksByAttributes(final String attributeNames[], final boolean[] ascending) {

        assert attributeNames.length == ascending.length;

//        for (TrackGroup tg : trackGroups) {
//            tg.sortByAttributes(attributeNames, ascending);
//        }
    }


    public void sortTracksByPosition(List<String> trackIds) {
//        for (TrackGroup tg : trackGroups) {
//            tg.sortByList(trackIds);
//        }
    }


    /**
     * Sort all groups (data and feature) by a computed score over a region.  The
     * sort is done twice (1) groups are sorted with the featureGroup, and (2) the
     * groups themselves are sorted.
     *
     * @param region
     * @param type
     */
    public void sortByRegionsScore(final RegionOfInterest region, final RegionScoreType type,
                                   final ReferenceFrame frame, List<String> sortedSamples) {

        sortGroupsByRegionScore(trackGroups, region, type, frame.getZoom(), frame.getName());
        //for (TrackGroup group : trackGroups) {
        // If there is a non-null linking attribute
        // Segregate tracks into 2 sub-groups, those matching the score type and those that do not
        //group.sortGroup(type, sortedSamples);
        //}
    }

    /**
     * Sort groups by a score (not the tracks within the group).
     *
     * @param groups
     * @param region
     * @param type
     * @param inzoom
     * @param frameName
     */
    private void sortGroupsByRegionScore(List<TrackGroup> groups,
                                         final RegionOfInterest region,
                                         final RegionScoreType type,
                                         int inzoom,
                                         final String frameName) {
        if ((groups != null) && (region != null) && !groups.isEmpty()) {
            final int zoom = Math.max(0, inzoom);
            final String chr = region.getChr();
            final int start = region.getStart();
            final int end = region.getEnd();
            Comparator<TrackGroup> c = (group1, group2) -> {
                float s1 = group1.getRegionScore(chr, start, end, zoom, type, frameName);
                float s2 = group2.getRegionScore(chr, start, end, zoom, type, frameName);
                // Use the Float comparator as it handles NaN.  Need to flip the order to make it descending
                return Float.compare(s2, s1);
            };
            Collections.sort(groups, c);
        }
    }

    public void removeTracks(Collection<? extends Track> tracksToRemove) {
        for (Track t : tracksToRemove) {
            if (t == this.track) {
                clearTracks();
                break;
            }
        }
    }

    /**
     * Remove, but do not dispose of, tracks.  Used by session reader
     */
    public void removeAllTracks() {
        track = null;
    }

    /**
     * Insert the selectedTracks collection either before or after the target and return true.
     *
     * @param selectedTracks
     * @param targetTrack
     * @param before
     */
    public void moveSelectedTracksTo(Collection<? extends Track> selectedTracks,
                                     Track targetTrack,
                                     boolean before) {

        if (selectedTracks.isEmpty()) {
            return;
        }

//        for (TrackGroup tg : trackGroups) {
//            if (tg.moveSelectedTracksTo(selectedTracks, targetTrack, before)) {
//                return;
//            }
//        }
    }


    public int getContentHeight() {
        return this.track.getContentHeight();
    }

    @Override
    public Dimension getPreferredSize() {
        return isVisible() ? new Dimension(100000, track.getContentHeight()) : new Dimension(0, 0);
    }

    @Override
    public Dimension getMaximumSize() {
        return getPreferredSize();
    }

    // Scrollable interface methods
    @Override
    public Dimension getPreferredScrollableViewportSize() {
        return getPreferredSize();
    }

    @Override
    public int getScrollableUnitIncrement(Rectangle visibleRect, int orientation, int direction) {
        return 16; // pixels per unit scroll
    }

    @Override
    public int getScrollableBlockIncrement(Rectangle visibleRect, int orientation, int direction) {
        if (orientation == SwingConstants.VERTICAL) {
            return visibleRect.height;
        }
        return visibleRect.width;
    }

    @Override
    public boolean getScrollableTracksViewportWidth() {
        // Return true to make the panel stretch to fit the viewport width
        return true;
    }

    @Override
    public boolean getScrollableTracksViewportHeight() {
        // Return false so the panel uses its preferred height and scrollbars appear when needed
        return false;
    }

    @Override
    public int getSnapshotHeight(boolean batch) {
        return getHeight();
    }

    public static TrackPanel getParentPanel(Track track) {
        for (TrackPanel panel : IGV.getInstance().getTrackPanels()) {
            for (TrackGroup group : panel.getGroups()) {
                if (group.contains(track)) {
                    return panel;
                }
            }
        }
        return null;
    }
}
