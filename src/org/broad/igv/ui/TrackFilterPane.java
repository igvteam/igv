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

package org.broad.igv.ui;

import org.broad.igv.track.Track;
import org.broad.igv.util.FilterElement;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * @author eflakes
 */
public class TrackFilterPane extends javax.swing.JPanel {

    private final int vgap;
    private boolean matchAll = true;
    private TrackFilter filter;
    protected List<String> itemList;
    protected String itemListLabel;
    private Component[] filterComponentBackup;
    JScrollPane scrollPane;

    public TrackFilterPane(List<String> items, String itemListLabel, TrackFilter filter) {

        setBackground(new java.awt.Color(255, 255, 255));
        setMinimumSize(new java.awt.Dimension(500, 100));
        vgap = 5;
        setLayout(new FlowLayout(FlowLayout.LEFT, 0, vgap));

        itemList = new ArrayList<String>(items);
        this.itemListLabel = itemListLabel;

        // Use the existing Filter if we have one otherwise create a new one
        if (filter != null) {

            this.filter = filter;
            FilterComponent filterComponent = null;

            Iterator iterator = filter.getFilterElements();
            while (iterator.hasNext()) {
                FilterElement element = (FilterElement) iterator.next();
                filterComponent = createFilterComponent(this, itemListLabel, itemList, element);
                filterComponent.displayMoreButton(false);
                add(filterComponent);
            }

            if (filterComponent != null) {
                filterComponent.displayMoreButton(true);
            }
        } else {
            this.filter = createNewFilter();
        }

        // Backup the initial state for restores
        backup();
        if (this.filter.isEmpty()) {
            more();
        }
    }

    public void setScrollPane(JScrollPane sp) {
        this.scrollPane = sp;
    }

    @Override
    public void setBounds(int x, int y, int width, int height) {
        System.out.println("Set bounds " + height);
        super.setBounds(x, y, width, height);    //To change body of overridden methods use File | Settings | File Templates.
    }

    @Override
    public void setBounds(Rectangle r) {
        System.out.println("Set bounds " + r.height);
        super.setBounds(r);    //To change body of overridden methods use File | Settings | File Templates.
    }

    @Override
    public void setPreferredSize(Dimension preferredSize) {
        System.out.println("Set pref size " + preferredSize.getHeight());
        super.setPreferredSize(preferredSize);    //To change body of overridden methods use File | Settings | File Templates.
    }

    public TrackFilter createNewFilter() {
        String name = "";
        return new TrackFilter(name, null);
    }

    public TrackFilter getFilter() {
        return filter;
    }

    public void clearTracks() {
        getFilter().clearTracks();
    }

    public void addTracks(List<Track> tracks) {
        getFilter().addTracks(tracks);
    }

    public FilterComponent createFilterComponent(TrackFilterPane filterPane,
                                                 String itemListLabel,
                                                 List<String> itemList,
                                                 FilterElement element) {

        return new TrackFilterComponent(filterPane, itemListLabel, itemList, element);
    }

    public void setMatchAll(boolean value) {
        matchAll = value;
    }

    public boolean getMatchAll() {
        return matchAll;
    }

    public void applyFilterMatching() {
        Component[] filterComponents = getComponents();
        for (Component filterComponent : filterComponents) {
            ((TrackFilterComponent) filterComponent).setMatchAll(matchAll);
        }
    }


    public void setItems(List<String> items) {

        if (items == null || items.size() < 1) {
            itemList = new ArrayList<String>();
            return;
        }
        itemList = new ArrayList<String>(items);
    }


    /**
     * Add another filter
     */
    public boolean more() {

        Window window = SwingUtilities.getWindowAncestor(this);

        if (!isFilterValid()) {
            JOptionPane.showMessageDialog(window,
                    "Some of the filter values are missing." +
                            "\nPlease enter all value before pressing ok.");
            return false;
        }

        // Add an unconfigured FilterComponent to the end of the list of
        // visual FilterComponents.
        if (itemListLabel == null) {
            itemListLabel = "When";
        }
        add(new TrackFilterComponent(this, itemListLabel, itemList, null));

        // Update prefered size
        int h = 0;
        for(Component c : getComponents()) {
            h += c.getHeight() + 2*vgap;
        }
        setPreferredSize(new Dimension(getWidth(), h));
        if(scrollPane != null) scrollPane.invalidate();

        return true;
    }

    public void save() {

        Component[] children = getComponents();
        for (Component child : children) {
            if (child instanceof FilterComponent) {
                ((FilterComponent) child).save();
            }
        }
    }

    public boolean isFilterValid() {

        boolean isValid = true;

        // Loop over each UI component that represent a single FilterElement
        // and check it's state for correctness
        Component[] children = getComponents();
        for (Component child : children) {
            if (child instanceof FilterComponent) {

                FilterComponent filterComponent = ((FilterComponent) child);

                String item = filterComponent.getItem();
                if (item == null || item.trim().equals(""))
                    return false;
                String operator = filterComponent.getComparisonOperator();
                if (operator == null || operator.trim().equals(""))
                    return false;
                String expectedValue = filterComponent.getExpectedValue();
                //if(expectedValue == null || expectedValue.trim().equals(""))
                //    return false;
            }
        }

        return isValid;
    }

    public void backup() {
        filterComponentBackup = getComponents();
    }

    public void adjustMoreAndBooleanButtonVisibility() {

        if (filterComponentBackup != null) {

            // Adjust the visibility of the More and the AND/OR Button
            Component[] components = this.getComponents();
            for (Component component : components) {
                ((FilterComponent) component).displayMoreButton(false);
            }

            int last = components.length - 1;
            ((FilterComponent) components[last]).displayMoreButton(true);
        }
    }

    public void restore() {

        if (filterComponentBackup != null) {

            // Remove all visual FilterComponents
            removeAll();

            // Remove all non-visual FilterElements
            filter.removeAll();

            // Adjust the visibility of the More and the AND/OR Button
            Component[] components = filterComponentBackup;
            for (Component component : components) {
                ((FilterComponent) component).displayMoreButton(false);

                // Put back the original FilterComponents
                add(component);

                // Put back the original non-visual FilterElements
                filter.add(((FilterComponent) component).getFilterElement());
            }

            if (components.length > 0) {

                // Adjust the visibility of the More
                int last = components.length - 1;
                ((FilterComponent) components[last]).displayMoreButton(true);
            }
            filterComponentBackup = null;
        }
    }


}
