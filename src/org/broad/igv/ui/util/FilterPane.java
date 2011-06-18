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

package org.broad.igv.ui.util;


import org.broad.igv.util.Filter;
import org.broad.igv.util.FilterElement;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


/**
 * @author eflakes
 */
abstract public class FilterPane extends javax.swing.JPanel {

    private Filter filter;
    protected List<String> itemList;
    protected String itemListLabel;
    private Component[] filterComponentBackup;


    /**
     * Creates new form FilterPane
     */
    public FilterPane(List<String> items, String itemListLabel, Filter filter) {

        initComponents();

        //setLayout( new SpringLayout() );
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

    abstract protected Filter createNewFilter();

    abstract protected FilterComponent createFilterComponent(
            FilterPane filterPane, String itemListLabel, List<String> items,
            FilterElement element);


    public Filter getFilter() {
        return filter;
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
        addDefaultFilterComponent();

        // Squeeze window to fit the new components
        if (window != null) {
            window.pack();
        }

        return true;
    }

    /**
     * Only for use with more() method.
     */
    abstract protected void addDefaultFilterComponent();

    @Override
    public Component add(Component component, int arg1) {
        if (component instanceof JButton) {

            if (((JButton) component).getText().equals("More...")) {
                return super.add(component, arg1);
            }
            return component;
        } else if (!(component instanceof FilterComponent))
            throw new RuntimeException("Invalid FilterComponent!");

        return super.add(component, arg1);
    }

    @Override
    public void add(Component component, Object arg1, int arg2) {

        if (component instanceof JButton) {

            if (((JButton) component).getText().equals("More...")) {
                super.add(component, arg1, arg2);
            }
        } else if (!(component instanceof FilterComponent))
            throw new RuntimeException("Invalid FilterComponent!");

        super.add(component, arg1, arg2);
    }

    @Override
    public void add(Component component, Object arg1) {

        if (component instanceof JButton) {

            if (((JButton) component).getText().equals("More...")) {
                super.add(component, arg1);
            }
        } else if (!(component instanceof FilterComponent))
            throw new RuntimeException("Invalid FilterComponent!");

        super.add(component, arg1);
    }

    @Override
    public Component add(Component component) {

        if (component instanceof JButton) {

            if (((JButton) component).getText().equals("More...")) {
                return super.add(component);
            }
            return component;
        } else if (!(component instanceof FilterComponent))
            throw new RuntimeException("Invalid FilterComponent!");

        return super.add(component);
    }

    @Override
    public Component add(String arg0, Component component) {

        if (component instanceof JButton) {

            if (((JButton) component).getText().equals("More...")) {
                return super.add(arg0, component);
            }
            return component;
        } else if (!(component instanceof FilterComponent))
            throw new RuntimeException("Invalid FilterComponent!");

        return super.add(arg0, component);

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

    protected void adjustMoreAndBooleanButtonVisibility() {

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


    /**
     * This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        setBackground(new java.awt.Color(255, 255, 255));
        setMinimumSize(new java.awt.Dimension(500, 100));
        setPreferredSize(new java.awt.Dimension(700, 100));
        setLayout(new java.awt.FlowLayout(java.awt.FlowLayout.LEFT, 5, 0));
    }// </editor-fold>//GEN-END:initComponents


    // Variables declaration - do not modify//GEN-BEGIN:variables
    // End of variables declaration//GEN-END:variables

}
