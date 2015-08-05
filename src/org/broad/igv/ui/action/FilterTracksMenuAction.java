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

import org.broad.igv.track.AttributeManager;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.TrackFilter;
import org.broad.igv.ui.TrackFilterPane;
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

        Integer response = showFilterTrackDialog(mainFrame.getMainFrame(), trackFilterPane, "Filter Tracks");

        if (response == null) {
            return;
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
            mainFrame.doRefresh();
        }

    }

    private Integer showFilterTrackDialog(Frame parent,
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
        dialog.setResizable(true);
        dialog.setBackground(Color.WHITE);
        dialog.getContentPane().setBackground(Color.WHITE);

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
                track.setVisible(true);
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
