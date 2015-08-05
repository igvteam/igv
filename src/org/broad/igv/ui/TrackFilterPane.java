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
