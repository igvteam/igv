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

import org.broad.igv.track.AttributeManager;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.TrackFilter;
import org.broad.igv.ui.TrackFilterPane;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.Filter;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.List;

/**
 * @author jrobinso
 */
public class FilterTracksMenuAction extends MenuAction {

    //static Logger log = Logger.getLogger(FilterTracksMenuAction.class);
    IGV mainFrame;

    private JCheckBox showAllTracksFilterCheckBox = new JCheckBox();

    private JCheckBox matchAllCheckBox = new JCheckBox();

    private JCheckBox matchAnyCheckBox = new JCheckBox();

    private TrackFilterPane trackFilterPane;

    public FilterTracksMenuAction(String label, int mnemonic, IGV mainFrame) {
        super(label, null, mnemonic);
        this.mainFrame = mainFrame;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        doFilterTracks();
    }

    private void doFilterTracks() {

        boolean previousDisableFilterState = showAllTracksFilterCheckBox.isSelected();

        boolean previousMatchAllState = matchAllCheckBox.isSelected();

        List<String> uniqueAttributeKeys = AttributeManager.getInstance().getAttributeNames();

        // Sort the attribute keys if we have any
        if (uniqueAttributeKeys != null) {
            //Collections.sort(uniqueAttributeKeys, AttributeManager.getInstance().getAttributeComparator());
        } else // If we have no attribute we can't display the
            // track filter dialog so say so and return
            if (uniqueAttributeKeys == null || uniqueAttributeKeys.isEmpty()) {

                MessageUtils.showMessage("No attributes found to use in a filter");
                return;
            }

        if (trackFilterPane == null) {
            trackFilterPane = new TrackFilterPane(uniqueAttributeKeys, "Show tracks whose attribute",
                    mainFrame.getSession().getFilter());

        } else {

            trackFilterPane.setItems(uniqueAttributeKeys);

            // Backup the initial state for restores
            trackFilterPane.backup();
            Filter filter = trackFilterPane.getFilter();
            if (filter == null || filter.isEmpty()) {
                trackFilterPane.more();
            }
        }

        trackFilterPane.clearTracks();
        trackFilterPane.addTracks(IGV.getInstance().getAllTracks());

        while (true) {

            Integer response = createFilterTrackDialog(mainFrame.getMainFrame(), trackFilterPane, "Filter Tracks");

            if (response == null) {
                continue;
            }

            if (response.intValue() == JOptionPane.CANCEL_OPTION) {

                // Restore previous filter state
                boolean disableFilterState = showAllTracksFilterCheckBox.isSelected();
                if (disableFilterState != previousDisableFilterState) {
                    showAllTracksFilterCheckBox.setSelected(previousDisableFilterState);
                }

                // Restore previous boolean match state
                boolean matchAllState = matchAllCheckBox.isSelected();
                if (matchAllState != previousMatchAllState) {
                    matchAllCheckBox.setSelected(previousMatchAllState);
                    matchAnyCheckBox.setSelected(!previousMatchAllState);
                }
                // Reset state
                trackFilterPane.restore();
                return;
            } else if (response.intValue() == JOptionPane.OK_OPTION) {

                filterTracks(trackFilterPane);
                break;
            }
        }

        // Update the state of the current tracks
        mainFrame.doRefresh();

    }

    private Integer createFilterTrackDialog(Frame parent,
                                            final TrackFilterPane trackFilterPane,
                                            String title) {

        JScrollPane scrollPane = new JScrollPane(trackFilterPane);
        scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
        scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

        int optionType = JOptionPane.OK_CANCEL_OPTION;
        int messageType = JOptionPane.PLAIN_MESSAGE;

        JPanel panel = new JPanel();
        panel.setLayout(new BorderLayout());

        JPanel filterHeaderPanel = new JPanel();
        filterHeaderPanel.setBackground(Color.WHITE);
        filterHeaderPanel.setLayout(new GridLayout(0, 1));
        filterHeaderPanel.add(new JLabel("For attributes that:"));

        ButtonGroup booleanButtonGroup = new ButtonGroup();
        booleanButtonGroup.add(matchAllCheckBox);
        booleanButtonGroup.add(matchAnyCheckBox);

        showAllTracksFilterCheckBox.setText("Show All Tracks");
        matchAllCheckBox.setText("Match all of the following");
        matchAnyCheckBox.setText("Match any of the following");
        boolean matchAll = trackFilterPane.getMatchAll();
        if (matchAll) {
            matchAllCheckBox.setSelected(true);
        } else {
            matchAnyCheckBox.setSelected(true);
        }

        matchAllCheckBox.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(java.awt.event.ActionEvent evt) {
                trackFilterPane.setMatchAll(true);
            }
        });
        matchAnyCheckBox.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(java.awt.event.ActionEvent evt) {
                trackFilterPane.setMatchAll(false);
            }
        });

        JPanel controls = new JPanel();
        FlowLayout layoutManager = new FlowLayout();
        layoutManager.setAlignment(FlowLayout.LEFT);
        controls.setLayout(layoutManager);
        controls.add(matchAllCheckBox);
        controls.add(matchAnyCheckBox);
        controls.add(showAllTracksFilterCheckBox);
        controls.setBackground(Color.WHITE);
        controls.setOpaque(true);
        filterHeaderPanel.add(controls);

        panel.setOpaque(true);
        panel.add(filterHeaderPanel, BorderLayout.NORTH);
        panel.add(scrollPane, BorderLayout.CENTER);

        final JOptionPane optionPane = new JOptionPane(panel, messageType, optionType);
        optionPane.setPreferredSize(new Dimension(700, 500));
        optionPane.setOpaque(true);
        optionPane.setBackground(Color.WHITE);
        optionPane.addPropertyChangeListener(
                JOptionPane.VALUE_PROPERTY,
                new PropertyChangeListener() {

                    public void propertyChange(PropertyChangeEvent e) {

                        Object value = e.getNewValue();
                        if (value instanceof Integer) {

                            int option = (Integer) value;
                            if (option == JOptionPane.OK_OPTION) {

                                if (trackFilterPane.isFilterValid()) {
                                    trackFilterPane.applyFilterMatching();
                                    trackFilterPane.save();
                                }
                            }
                        }
                    }
                });

        JDialog dialog = optionPane.createDialog(parent, title);
        dialog.setBackground(Color.WHITE);
        dialog.getContentPane().setBackground(Color.WHITE);

        Dimension maximumSize = new Dimension(dialog.getSize().width, 100);
        if (maximumSize != null) {
            dialog.setMaximumSize(maximumSize);
        }

        Component[] children = optionPane.getComponents();
        if (children != null) {
            for (Component child : children) {
                child.setBackground(Color.WHITE);
            }

        }

        dialog.pack();
        dialog.setVisible(true);

        Object selectedValue = optionPane.getValue();
        if (selectedValue == null) {
            return JOptionPane.CANCEL_OPTION;
        } else if (((Integer) selectedValue).intValue() == JOptionPane.OK_OPTION) {
            if (!trackFilterPane.isFilterValid() &&
                    !showAllTracksFilterCheckBox.isSelected()) {
                JOptionPane.showMessageDialog(
                        parent,
                        "Some of the filter values are missing." +
                                "\nPlease enter all value before pressing ok.");

                selectedValue = null;
            }
        }
        return ((Integer) selectedValue);
    }

    private void filterTracks(TrackFilterPane trackFilterPane) {

        boolean showAllTracks = showAllTracksFilterCheckBox.isSelected();
        if (showAllTracks) {

            List<Track> tracks = IGV.getInstance().getAllTracks();
            for (Track track : tracks) {
                track.setVisible(showAllTracks);
            }

        } else {
            TrackFilter filter = trackFilterPane.getFilter();
            IGV.getInstance().getSession().setFilter(filter);
            // Evaluate the filter elements
            filter.evaluate();
        }

    }

    public void resetTrackFilter() {
        trackFilterPane = null;
        IGV.getInstance().getSession().setFilter(null);
        setFilterShowAllTracks(false);
    }

    public void setFilterShowAllTracks(boolean value) {
        if (showAllTracksFilterCheckBox != null) {
            showAllTracksFilterCheckBox.setSelected(value);
        }
    }

    public JCheckBox getShowAllTracksFilterCheckBox() {
        return showAllTracksFilterCheckBox;
    }

    public void updateTrackFilter() {

        TrackFilter trackFilter = IGV.getInstance().getSession().getFilter();

        if (trackFilter == null) {
            return;
        }

        List<String> uniqueAttributeKeys = AttributeManager.getInstance().getAttributeNames();


        // If we have no attribute we can't display the
        // track filter dialog so say so and return
        if (uniqueAttributeKeys == null || uniqueAttributeKeys.isEmpty()) {

            MessageUtils.showMessage("No attributes found to use in a filter");
            return;
        }

        trackFilterPane = new TrackFilterPane(uniqueAttributeKeys, "Show tracks whose attribute", trackFilter);
        trackFilterPane.clearTracks();
        trackFilterPane.addTracks(IGV.getInstance().getAllTracks());

        // Evaluate the filter elements
        trackFilter.evaluate();
    }

    public void setFilterMatchAll(boolean value) {
        if (trackFilterPane != null) {
            trackFilterPane.setMatchAll(value);
        }
    }

    public boolean isFilterMatchAll() {
        if (trackFilterPane != null) {
            return trackFilterPane.getMatchAll();
        }
        return false;
    }
}
