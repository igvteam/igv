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
import java.util.stream.Collectors;

/**
 * @author eflakes
 */
public class TrackPanel extends JPanel implements Paintable {

    MainPanel mainPanel;

    private static Logger log = LogManager.getLogger(TrackPanel.class);

    private String name = null;
    private TrackNamePanel namePanel;
    private AttributePanel attributePanel;
    private DataPanelContainer dataPanelContainer;
    private List<Track> tracks;
    transient int lastHeight = 0;

    transient int visibleHeight;
    private TrackPanelScrollPane scrollPane;

    /**
     * Constructs ...
     *
     * @param name
     */
    public TrackPanel(String name, MainPanel mainPanel) {
        setLayout(null); //new TrackPanelLayout());
        this.mainPanel = mainPanel;
        this.name = name;
        tracks = Collections.synchronizedList(new LinkedList<Track>());
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
        return scrollPane;
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

    /**
     * Override setVisible to also hide the containing scroll pane.
     *
     * @param aFlag true to make the component visible; false to
     *              make it invisible
     */
    @Override
    public void setVisible(boolean aFlag) {
        //super.setVisible(aFlag);
        getScrollPane().setVisible(aFlag);
        IGV.getInstance().getMainPanel().revalidateTrackPanels();
    }

    public void paintOffscreen(Graphics2D g, Rectangle rect, boolean batch) {

        g.setColor(Color.black);

        Component[] children = getComponents();

        for (Component component : children) {
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
    public boolean hasTracks() {
        return !tracks.isEmpty();
    }

    public int getVisibleTrackCount() {
        int count = 0;
        for (Track t : tracks) {
            if (t != null && t.isVisible()) {
                count++;
            }
        }
        return count;
    }

    public List<Track> getTracks() {
        return tracks;
    }

    public List<Track> getVisibleTracks() {
        return tracks.stream().filter(t -> t.isVisible()).collect(Collectors.toList());
    }

    public void clearTracks() {
        tracks.clear();
    }


    @Deprecated
    public boolean fitTracksToPanel() {
        return true;
    }

    public void addTrack(Track track) {
        tracks.add(track);
    }

    public void addTracks(Collection<? extends Track> tracks) {
        for (Track t : tracks) {
            addTrack(t);
        }
    }

    public void reset() {
        clearTracks();
    }


    public void removeTracks(Collection<? extends Track> tracksToRemove) {
        tracks.removeAll(tracksToRemove);
    }

    /**
     * Remove, but do not dispose of, tracks.  Used by session reader
     */
    public void removeAllTracks() {
        tracks.clear();
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

        int index = (targetTrack == null ? tracks.size() : tracks.indexOf(targetTrack));
        if (index > 0) {

            if (!before) {
                index = index + 1;
            }

            // 1. Divdide the target list up into 2 parts, one before the index and one after
            List<Track> beforeList = new ArrayList(tracks.subList(0, index));
            List<Track> afterList = new ArrayList(tracks.subList(index, tracks.size()));

            // 2.  Remove the selected tracks from anywhere they occur
            beforeList.removeAll(selectedTracks);
            afterList.removeAll(selectedTracks);

            // 3. Now insert the selected tracks
            tracks.clear();
            tracks.addAll(beforeList);
            tracks.addAll(selectedTracks);
            tracks.addAll(afterList);
        }

    }

    /**
     * Prefered size == size to view all content without scrolling.
     *
     * @return
     */
    @Override
    public Dimension getPreferredSize() {
        Dimension dim = super.getPreferredSize();
        dim.height = getPreferredPanelHeight();
        return dim;
    }

//    public void updatePreferredSize() {
//
//        int height = 0;
//        for (Track t : tracks) {
//            height += t.getContentHeight();
//        }
//        TrackPanelScrollPane sp = getScrollPane();
//        Insets insets2 = sp.getInsets();
//        int dy = insets2.top + insets2.bottom;
//        sp.setPreferredSize(new Dimension(1000, height + dy));
//        mainPanel.revalidate();
//    }

    public int getPreferredPanelHeight() {

        int h = 0;
        for (Track track : tracks) {
            if (track != null && track.isVisible()) {
                h += track.getContentHeight();
            }
        }
        return h;
    }


    /**
     * Return the sum of track.height values for all tracks associated with this panel.  Normally this is a single
     * track.  The "track.height" property corresponds to the height from the users perspective, i.e. the height
     * of the associated scroll pane
     */
    public int getTotalTrackHeight() {
        return tracks.stream().collect(Collectors.summingInt(t -> t.isVisible() ? t.getHeight() : 0));
    }

    public boolean isHeightChanged() {
        int height = getPreferredPanelHeight();
        boolean change = height != lastHeight;
        lastHeight = height;
        return change;
    }


    @Override
    public Dimension getMinimumSize() {
        return new Dimension(0, 50);
    }

    public Dimension getMaximumSize() {
        return new Dimension(32767, 500);  // 32767 is hardcoded in Swing
    }


    public void sortTracksByAttributes(final String attributeNames[], final boolean[] ascending) {
        assert attributeNames.length == ascending.length;
        tracks = TrackPanelHelper.sortByAttributes(tracks, attributeNames, ascending);
    }

    public void sortTracksByPosition(List<String> trackIds) {
        tracks = TrackPanelHelper.sortByList(tracks, trackIds);
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

//        sortGroupsByRegionScore(trackGroups, region, type, frame.getZoom(), frame.getName());
//        for (TrackGroup group : trackGroups) {
//            // If there is a non-null linking attribute
//            // Segregate tracks into 2 sub-groups, those matching the score type and those that do not
//            group.sortGroup(type, sortedSamples);
//        }
//
        tracks = TrackPanelHelper.sortByRegionScore(tracks, type, sortedSamples);
    }


    public void setScrollPane(TrackPanelScrollPane sp) {
        this.scrollPane = sp;
    }

}
