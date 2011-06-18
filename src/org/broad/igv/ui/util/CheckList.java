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

package org.broad.igv.ui.util;

import org.broad.igv.util.Utilities;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.*;
import java.util.List;

/**
 * @author eflakes
 */
public class CheckList extends JPanel {

    private boolean editing = false;
    private boolean isDirty = false;
    private HashSet<String> itemSet = new HashSet<String>();
    private boolean defaultState = true;
    private JLabel header = null;
    private JPanel itemPanel = new JPanel();
    private GridLayout layout;

    // Backup of state
    private Component[] componentState = new Component[0];

    public CheckList(boolean defaultState) {
        this(defaultState, null);
    }

    public CheckList(boolean defaultState, String headerText) {

        setLayout(new BorderLayout());
        this.defaultState = defaultState;
        layout = new GridLayout(0, 1);
        itemPanel.setLayout(layout);

        if (headerText != null) {
            header = new JLabel(headerText);
            add(header, BorderLayout.NORTH);
            add(itemPanel, BorderLayout.CENTER);
        } else {
            add(itemPanel);
        }

        JPanel checkClearPanel = new JPanel();
        checkClearPanel.setLayout(new FlowLayout());

        JButton checkAllButton = new JButton("Select All");
        checkAllButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent actionEvent) {
                setStateForAll(true);
            }
        });

        JButton clearAllButton = new JButton("Clear All");
        clearAllButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent actionEvent) {
                setStateForAll(false);
            }
        });


        checkClearPanel.add(checkAllButton);
        checkClearPanel.add(clearAllButton);
        add(checkClearPanel, BorderLayout.SOUTH);
    }


    private void setStateForAll(boolean selected) {
        for(Component c : itemPanel.getComponents()) {
            if(c instanceof JCheckBox) {
                ((JCheckBox) c).setSelected(selected);
            }
        }
        isDirty = true;
    }

    public void addItems(List<String> itemStrings) {

        addItems(itemStrings, defaultState);
    }

    public void addItems(List<String> itemStrings, boolean isChecked) {

        // Layout so that we have a maximum of 20 items per column
        int nColumns = (itemStrings.size() / 20) + 1;
        layout.setColumns(nColumns);

        for (String itemText : itemStrings) {
            addItem(itemText, isChecked);
        }
        JButton checkAllButton = new JButton("Check All");
        itemPanel.add(checkAllButton);
        isDirty = true;
    }

    public void addItem(String itemText, Boolean isChecked) {

        if (itemSet.contains(itemText.trim())) {
            return;
        }

        JCheckBox checkBox = new JCheckBox();
        if (itemText != null) {
            checkBox.setText(itemText);
        }
        if (isChecked == null) {
            checkBox.setSelected(defaultState);
        } else {
            checkBox.setSelected(isChecked);
        }
        itemPanel.add(checkBox);
        itemSet.add(itemText.trim());
        isDirty = true;
    }

    public void clear() {
        itemPanel.removeAll();
        itemSet.clear();
        isDirty = true;
    }

    /**
     * Make sure the current state of the checkboxes are used. Updating
     * is not allowed while in edit mode.
     */
    public void update() {

        // Can't restorePersistentState while in the middle of editing
        if (isEditing()) {
            return;
        }

        Component[] components = getItemComponents();
        int length = components.length;
        componentState = new Component[length];
        int i = 0;

        for (Component component : itemPanel.getComponents()) {

            boolean isChecked = ((JCheckBox) component).isSelected();
            String text = ((JCheckBox) component).getText();

            JCheckBox checkBox = new JCheckBox();
            checkBox.setSelected(isChecked);
            checkBox.setText(text);
            componentState[i++] = checkBox;
        }
    }

    public boolean isEditing() {
        return editing;
    }

    public void setEditing(boolean editing) {
        this.editing = editing;
    }

    public HashSet<String> getSelectedItems() {

        // If there has been a modification we need to restorePersistentState first
        if (isDirty) {
            update();
        }

        HashSet<String> set = new LinkedHashSet<String>();

        for (Component component : componentState) {

            if (component instanceof JCheckBox) {
                JCheckBox checkBox = ((JCheckBox) component);

                if (!checkBox.isSelected())
                    continue;

                String text = checkBox.getText().trim();
                if (text != null && !text.equals("")) {
                    set.add(text);
                }
            }
        }
        return set;
    }

    public HashSet<String> getUnselectedItems() {

        // If there has been a modification we need to restorePersistentState first
        if (isDirty) {
            update();
        }

        HashSet<String> set = new LinkedHashSet<String>();

        for (Component component : componentState) {

            if (component instanceof JCheckBox) {
                JCheckBox checkBox = ((JCheckBox) component);

                if (checkBox.isSelected())
                    continue;

                String text = checkBox.getText().trim();
                if (text != null && !text.equals("")) {
                    set.add(text);
                }
            }
        }
        return set;
    }

    public void deselectItems(Collection<String> strings) {

        for (Component component : componentState) {

            if (component instanceof JCheckBox) {
                JCheckBox checkBox = ((JCheckBox) component);
                if (checkBox.isSelected()) {
                    String text = checkBox.getText().trim();
                    if (strings.contains(text)) {
                        checkBox.setSelected(false);
                    }
                }
            }
        }
    }

    protected Component[] getItemComponents() {
        return itemPanel.getComponents();
    }

    public void sort() {

        List<String> temp = new ArrayList<String>();
        for (String item : itemSet) {
            temp.add(item.toLowerCase());
        }

        // Do the sort
        Collections.sort(temp, Utilities.getNumericStringComparator());

        Component[] checkBoxes = itemPanel.getComponents();
        itemPanel.removeAll();
        for (int i = 0; i < temp.size(); i++) {
            for (Component checkBox : checkBoxes) {
                String name = ((JCheckBox) checkBox).getText();
                if (name.equalsIgnoreCase(temp.get(i))) {
                    itemPanel.add(checkBox);
                    break;
                }
            }
        }
        isDirty = true;
    }

    public void cancelChanges() {

        // Can't cancel changes if not editing
        if (!isEditing()) {
            return;
        }

        // Remove the changes
        itemPanel.removeAll();

        // Put the previous state back
        for (Component item : componentState) {
            itemPanel.add(item);
        }
    }
}
