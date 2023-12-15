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

package org.broad.igv.ui.panel;


import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackGroup;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;

import javax.swing.*;
import java.awt.*;
import java.util.List;
import java.util.*;

/**
 * @author eflakes
 */
public class TrackPanel extends JPanel implements Paintable  { //} implements Scrollable {

    MainPanel mainPanel;

    private static Logger log = LogManager.getLogger(TrackPanel.class);

    private String name = null;
    private TrackNamePanel namePanel;
    private AttributePanel attributePanel;
    private DataPanelContainer dataPanelContainer;
    private String groupAttribute;
    private List<TrackGroup> trackGroups;
    transient int lastHeight = 0;

    /**
     * Constructs ...
     *
     * @param name
     */
    public TrackPanel(String name, MainPanel mainPanel) {
        setLayout(null); //new TrackPanelLayout());
        this.mainPanel = mainPanel;
        this.name = name;
        TrackGroup nullGroup = new TrackGroup();
        nullGroup.setDrawBorder(false);
        trackGroups = Collections.synchronizedList(new LinkedList<TrackGroup>());
        trackGroups.add(nullGroup);
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

    public int getViewportHeight() {

        Container parent = getParent();
        return parent == null ? 0 : parent.getHeight();
    }

    public TrackPanelScrollPane getScrollPane() {

        TrackPanelScrollPane scollpane = null;
        Container parent = getParent();

        if (parent instanceof JViewport) {
            scollpane = (TrackPanelScrollPane) ((JViewport) parent).getParent();
        }
        return scollpane;
    }

    @Override
    public void doLayout() {

        synchronized (getTreeLock()) {

            int h = getHeight(); //getPreferredSize().height;
            Component[] children = getComponents();

            int nw = mainPanel.getNamePanelWidth();

            Component namePanel = children[0];
            Component attributePanel = children[1];
            Component dataPanel = children[2];
            namePanel.setBounds(mainPanel.getNamePanelX(), 0, nw, h);

            attributePanel.setBounds(mainPanel.getAttributePanelX(), 0, mainPanel.getAttributePanelWidth(), h);  // Attributes
            dataPanel.setBounds(mainPanel.getDataPanelX(), 0, mainPanel.getDataPanelWidth(), h);

            dataPanel.doLayout();
        }
    }

    public void paintOffscreen(Graphics2D g, Rectangle rect, boolean batch) {

        g.setColor(Color.black);

        Component[] children = getComponents();

        for(Component component : children) {
            if (component instanceof Paintable) {
                if (component.getWidth() > 0) {
                    Rectangle clipRect = new Rectangle(0, rect.y, component.getWidth(), rect.height);
                    Graphics2D graphics = (Graphics2D) g.create();
                    graphics.translate(component.getX(), component.getY());
                    //  graphics.setClip(clipRect);
                    ((Paintable) component).paintOffscreen(graphics, clipRect, batch);
                    graphics.dispose();
                }
            }
        }
    }

    @Override
    public int getSnapshotHeight(boolean batch) {
        return getHeight();
    }

    public void createDataPanels() {
        dataPanelContainer.createDataPanels();
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
        for (TrackGroup tg : trackGroups) {
            if (tg.getVisibleTracks().size() > 0) {
                return true;
            }
        }
        return false;
    }

    public int getVisibleTrackCount() {
        int count = 0;
        for (TrackGroup tg : trackGroups) {
            for (Track t : tg.getVisibleTracks()) {
                if (t != null && t.isVisible()) {
                    count++;
                }
            }
        }
        return count;
    }

    public boolean isHeightChanged() {
        int height = getPreferredPanelHeight();
        boolean change = height != lastHeight;
        lastHeight = height;
        return change;
    }

    public List<Track> getTracks() {
        ArrayList<Track> tracks = new ArrayList();
        for (TrackGroup tg : trackGroups) {
            tracks.addAll(tg.getTracks());
        }
        return tracks;
    }

    public void clearTracks() {

        for (Track t : getTracks()) {
            t.unload();
        }
        groupAttribute = null;
        trackGroups.clear();
    }


    public boolean fitTracksToPanel() {
        DataPanelContainer dataPanel = this.getScrollPane().getDataPanel();
        boolean success = true;

        int availableHeight = dataPanel.getVisibleHeight();
        int visibleTrackCount = 0;

        // Process data tracks first
        Collection<TrackGroup> groups = dataPanel.getTrackGroups();


        // Count visible tracks.
        for (TrackGroup group : groups) {
            List<Track> tracks = group.getVisibleTracks();
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
                List<Track> tracks = group.getVisibleTracks();
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

    /**
     * Add a track to this panel.  If tracks are grouped, search for correct group, or make a new one if not found.
     *
     * @param track
     */
    public void addTrack(Track track) {
        log.debug("Adding track " + track.getName() + " to panel " + getName());
        String groupName = (groupAttribute == null ? null : track.getAttributeValue(groupAttribute));
        boolean foundGroup = false;
        for (TrackGroup tg : trackGroups) {
            if (groupAttribute == null || groupName == null || tg.getName().equals(groupName)) {
                tg.add(track);
                foundGroup = true;
                break;
            }
        }
        if (!foundGroup) {
            TrackGroup newGroup = new TrackGroup(groupName);
            newGroup.add(track);
            if (groupAttribute == null) {
                newGroup.setDrawBorder(false);
            }
            trackGroups.add(newGroup);
        }
    }

    public void addTracks(Collection<? extends Track> tracks) {
        for (Track t : tracks) {
            addTrack(t);
        }
    }

    public void moveGroup(TrackGroup group, int index) {

        if (index > trackGroups.indexOf(group)) {
            index--;
        }
        trackGroups.remove(group);
        if (index >= trackGroups.size()) {
            trackGroups.add(group);
        } else {
            trackGroups.add(index, group);
        }
    }


    public void reset() {
        this.groupAttribute = null;
        clearTracks();
    }

    /**
     * Rebuild group list for supplied attribute.
     *
     * @param attribute
     */
    public void groupTracksByAttribute(String attribute) {

        this.groupAttribute = attribute;
        List<Track> tracks = getTracks();
        trackGroups.clear();

        if (attribute == null || attribute.length() == 0) {
            TrackGroup nullGroup = new TrackGroup();
            nullGroup.addAll(tracks);
            nullGroup.setDrawBorder(false);
            trackGroups.add(nullGroup);
        } else {
            Map<String, TrackGroup> groupMap = new HashMap();
            for (Track track : tracks) {
                String attributeValue = track.getAttributeValue(attribute);

                if (attributeValue == null) {
                    attributeValue = "";
                }

                TrackGroup group = groupMap.get(attributeValue);

                if (group == null) {
                    group = new TrackGroup(attributeValue);
                    groupMap.put(attributeValue, group);
                    trackGroups.add(group);
                }
                group.add(track);
            }
        }
    }

    public void sortTracksByAttributes(final String attributeNames[], final boolean[] ascending) {

        assert attributeNames.length == ascending.length;

        for (TrackGroup tg : trackGroups) {
            tg.sortByAttributes(attributeNames, ascending);
        }
    }


    public void sortTracksByPosition(List<String> trackIds) {
        for (TrackGroup tg : trackGroups) {
            tg.sortByList(trackIds);
        }

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
        for (TrackGroup group : trackGroups) {
            // If there is a non-null linking attribute
            // Segregate tracks into 2 sub-groups, those matching the score type and those that do not
            group.sortGroup(type, sortedSamples);
        }
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
        for (TrackGroup tg : trackGroups) {
            tg.removeTracks(tracksToRemove);
        }
    }

    /**
     * Remove, but do not dispose of, tracks.  Used by session reader
     */
    public void removeAllTracks() {
        trackGroups.clear();
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

        for (TrackGroup tg : trackGroups) {
            if (tg.moveSelectedTracksTo(selectedTracks, targetTrack, before)) {
                return;
            }
        }
    }

    public int getPreferredPanelHeight() {

        int height = 0;

        Collection<TrackGroup> groups = getGroups();

        if (groups.size() > 1) {
            height += UIConstants.groupGap;
        }

        synchronized (groups) {

            for (Iterator<TrackGroup> groupIter = groups.iterator(); groupIter.hasNext(); ) {
                TrackGroup group = groupIter.next();
                if (group != null && group.isVisible()) {
                    if (groups.size() > 1) {
                        height += UIConstants.groupGap;
                    }
                    height += group.getHeight();
                }
            }
        }

        return height;
    }

    /**
     * Prefered size == size to view all content without scrolling.
     * @return
     */
    @Override
    public Dimension getPreferredSize() {
        Dimension dim = super.getPreferredSize();
        dim.height = getPreferredPanelHeight();
        return dim;
    }

    public void updatePreferredSize() {

        int height = 0;
        int minimumHeight = 0;
        for(TrackGroup trackGroup : this.trackGroups) {
            for(Track t : trackGroup.getVisibleTracks());
            height += trackGroup.getHeight();
                    }

        TrackPanelScrollPane sp = getScrollPane();
        Insets insets2 = sp.getInsets();
        int dy = insets2.top + insets2.bottom;
        sp.setPreferredSize(new Dimension(1000, height + dy));
        mainPanel.revalidate();
    }


    @Override
    public Dimension getMinimumSize() {
        return new Dimension(0, 50);
    }

    public Dimension getMaximumSize() {
        return new Dimension(32767, 500);  // 32767 is hardcoded in Swing
    }

    public void addTrackGroup(TrackGroup trackGroup) {
        trackGroups.add(trackGroup);
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
