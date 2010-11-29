/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.WaitCursorManager;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 * @author eflakes
 */
public class PropertyDialog extends OkCancelDialog {

    // General

    /**
     * Field description
     */
    final static public String CLICK_ITEM_TO_EDIT_TOOLTIP = "Click this item bring up its editor";

    static public enum PreferenceType {
        TEXT, BOOLEAN, COLOR
    }

    /**
     * Class description
     *
     * @author Enter your name here...
     * @version Enter version here..., 08/11/24
     */
    static public class PreferenceDescriptor {

        private PreferenceType valueType = PreferenceType.TEXT;
        private String valueKey;
        private Object defaultValue;

        /**
         * Constructs ...
         *
         * @param valueKey
         * @param valueType
         */
        public PreferenceDescriptor(String valueKey, PreferenceType valueType) {
            this(valueKey, valueType, null);
        }

        /**
         * Constructs ...
         *
         * @param valueKey
         * @param valueType
         * @param defaultValue
         */
        public PreferenceDescriptor(String valueKey, PreferenceType valueType,
                                    Object defaultValue) {

            this.valueKey = valueKey;
            this.valueType = valueType;
            this.defaultValue = defaultValue;
        }

        PreferenceType getType() {
            return valueType;
        }

        String getValueKey() {
            return valueKey;
        }
    }


    private PropertyManager propertyManager;
    private GroupContentPane content;

    /**
     * Creates new form PropertyDialog
     *
     * @param propertyManager
     * @param labelTextToKey
     * @param parent
     * @param modal
     */
    public PropertyDialog(PropertyManager propertyManager,
                          LinkedHashMap<String, PreferenceDescriptor> labelTextToKey,
                          Dialog parent, boolean modal) {

        super(parent, modal);
        initialize(propertyManager, labelTextToKey);
    }

    /**
     * Creates new form PropertyDialog
     *
     * @param propertyManager
     * @param labelTextToKey
     * @param parent
     * @param modal
     */
    public PropertyDialog(PropertyManager propertyManager,
                          LinkedHashMap<String, PreferenceDescriptor> labelTextToKey, Frame parent,
                          boolean modal) {

        super(parent, modal);
        initialize(propertyManager, labelTextToKey);
    }

    private void initialize(PropertyManager propertyManager,
                            LinkedHashMap<String, PreferenceDescriptor> labelTextToKey) {
        this.propertyManager = propertyManager;
        configureProperties(labelTextToKey);
        pack();
    }

    private void configureProperties(LinkedHashMap<String, PreferenceDescriptor> labelTextToKey) {

        // Get the area where our content can be displayed
        JPanel contentPane = getDialogPanel();
        contentPane.setLayout(new BorderLayout());

        content = getContent(labelTextToKey);

        JScrollPane scrollPane = new JScrollPane(content);
        scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
        scrollPane.setBorder(null);
        contentPane.add(scrollPane);
    }

    private GroupContentPane getContent(Map<String, PreferenceDescriptor> labelTextToKey) {

        GroupContentPane contentPanel = new GroupContentPane(new SpringLayout());

        // Property Group Settings
        int verticalGroupSpacing = 1;
        int groupHeight = 25;
        int groupWidth = 600;

        // Property Panel Settings
        int propertyPanelWidth = 600;
        int propertyPanelHeight = verticalGroupSpacing;

        // Layout configuration
        Spring topSpring = Spring.constant(verticalGroupSpacing);
        verticalGroupSpacing = groupHeight + 10;
        Spring increment = Spring.constant(verticalGroupSpacing);
        Spring centerSpring = Spring.constant((propertyPanelWidth - groupWidth) / 2);

        Set<Entry<String, PreferenceDescriptor>> set = labelTextToKey.entrySet();
        for (Entry<String, PreferenceDescriptor> entry : set) {

            String labelText = entry.getKey();
            PreferenceDescriptor preferenceDescriptor = entry.getValue();
            String externalPreferenceKey = preferenceDescriptor.getValueKey();
            PreferenceType preferenceType = preferenceDescriptor.getType();
            String defaultType = null;

            // If Preference type is boolean default to a boolean type
            if (preferenceType.equals(PreferenceType.BOOLEAN)) {
                defaultType = "false";
            }

            // VALUE LOOKUP: Use the value key to look up the value
            // from the preference store
            String valueText = get(externalPreferenceKey, defaultType);

            Dimension groupSize = new Dimension(groupWidth, groupHeight);

            if (valueText != null) {
                valueText = valueText.trim();
                if (valueText.equalsIgnoreCase("true") || valueText.equalsIgnoreCase("false")) {
                    preferenceType = PreferenceType.BOOLEAN;
                }
            }

            // Create the constraint
            SpringLayout.Constraints constraints = new SpringLayout.Constraints();
            constraints.setConstraint(SpringLayout.NORTH, topSpring);
            constraints.setConstraint(SpringLayout.WEST, centerSpring);

            // Create a group that contains the label and the value
            Group group = new Group(externalPreferenceKey, labelText, valueText, groupSize,
                    preferenceType, preferenceDescriptor.defaultValue);
            contentPanel.add(group, constraints);

            topSpring = Spring.sum(topSpring, increment);

            // calculate the need content panel height which depends
            // on the number of items in it
            propertyPanelHeight += verticalGroupSpacing;
        }

        Dimension contentSize = new Dimension(propertyPanelWidth, propertyPanelHeight + 10);
        contentPanel.setSize(contentSize);
        contentPanel.setPreferredSize(contentSize);
        contentPanel.setMinimumSize(contentSize);

        return contentPanel;
    }

    /**
     * Method description
     *
     * @param event
     * @return
     */
    public boolean okButtonClicked(java.awt.event.ActionEvent event) {

        // Read all items from the dialog then put them in the preference object
        Set<Entry<String, String>> set = content.getData().entrySet();
        for (Entry<String, String> entry : set) {

            String key = entry.getKey();
            String value = entry.getValue();

            if ((value == null) || value.trim().equals("")) {
                remove(key);
                continue;
            }
            put(key, value);
        }
        return true;
    }

    /**
     * Method description
     *
     * @param event
     * @return
     */
    public boolean cancelButtonClicked(java.awt.event.ActionEvent event) {
        return true;
    }

    protected PropertyManager getPreferences() {
        return propertyManager;
    }

    protected void put(String key, String value) {
        propertyManager.put(key, value);
    }

    protected String get(String key, String defaultValue) {
        return propertyManager.get(key, defaultValue);
    }

    protected void remove(String key) {
        propertyManager.remove(key);
    }

    /**
     * Panel represents a Colot type
     */
    static private class ColorPanel extends JPanel {

        private WaitCursorManager.CursorToken token;

        /**
         * Constructs ...
         */
        public ColorPanel() {

            addMouseListener(new MouseAdapter() {

                @Override
                public void mouseEntered(MouseEvent e) {

                    setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
                }

                @Override
                public void mouseExited(MouseEvent e) {
                }


                @Override
                public void mouseClicked(final MouseEvent e) {
                    final Color color = UIUtilities.showColorChooserDialog("", null);

                    if (color != null) {

                        JPanel panel = (JPanel) e.getSource();
                        panel.setBackground(color);
                        panel.repaint();
                    }
                }
            });

            UIUtilities.invokeOnEventThread(new Runnable() {
                public void run() {
                    ColorPanel.this.setToolTipText(UIConstants.CLICK_ITEM_TO_EDIT_TOOLTIP);
                }
            });
        }
    }


    /**
     * Handle getting the new values to be returned to called
     */
    static private class GroupContentPane extends JPanel {

        final static long serialVersionUID = 1198458423465L;

        /**
         * Constructs ...
         *
         * @param layoutManager
         */
        public GroupContentPane(LayoutManager layoutManager) {
            super(layoutManager);
        }

        /**
         * Get the new value set on the property dialog,
         *
         * @return
         */
        public Map<String, String> getData() {

            Map<String, String> keyAndValue = new LinkedHashMap<String, String>();

            // The children of this pane should be a bunch of Groups
            Component[] components = getComponents();
            for (Component component : components) {
                Group group = (Group) component;
                keyAndValue.put(group.getExternalPreferenceKey(), group.getValue());
            }
            return keyAndValue;
        }
    }


    /**
     * Group for label and value
     */
    static private class Group extends JPanel {

        final static long serialVersionUID = 1198458423465L;

        private String externalPreferenceKey;
        private JLabel label;
        private JComponent value;
        private StringBuffer scratchBuffer = new StringBuffer();

        /**
         * Constructs ...
         */
        public Group() {
            super();
            setLayout(new FlowLayout(FlowLayout.LEFT));
        }

        /**
         * Constructs ...
         *
         * @param preferenceKey
         * @param labelText
         * @param valueText
         * @param size
         * @param valueType
         * @param defaultValue
         */
        public Group(String preferenceKey, String labelText, String valueText, Dimension size,
                     PreferenceType valueType, Object defaultValue) {

            int labelWidth = (int) (size.width / 2.5);
            int labelHeight = size.height;
            int valueHeight = labelHeight;

            Dimension dimension = null;

            this.externalPreferenceKey = preferenceKey;

            // Label
            dimension = new Dimension(labelWidth, labelHeight);
            label = new JLabel(labelText);
            label.setHorizontalAlignment(SwingConstants.LEFT);
            label.setSize(dimension);
            label.setPreferredSize(dimension);
            label.setMinimumSize(dimension);
            label.setOpaque(true);
            add(label);


            if (valueType.equals(PreferenceType.BOOLEAN)) {

                valueText = valueText.toLowerCase();
                dimension = new Dimension(labelHeight, valueHeight);
                value = new JCheckBox();

                // If no value was set but we have a default use it
                if ((valueText == null) && (defaultValue != null)) {

                    ((JCheckBox) value).setSelected(
                            new Boolean(((String) defaultValue)).booleanValue());
                } else {

                    ((JCheckBox) value).setSelected(new Boolean(valueText).booleanValue());
                }

                value.setSize(dimension);
                value.setPreferredSize(dimension);
                value.setMinimumSize(dimension);
                add(value);
            } else if (valueType.equals(PreferenceType.TEXT)) {

                int valueWidth = (int) (labelWidth * 1.4);
                dimension = new Dimension(valueWidth, valueHeight);

                // If no value was set but we have a default use it
                if ((valueText == null) && (defaultValue != null)) {

                    value = new JTextField((String) defaultValue);
                } else {

                    value = new JTextField(valueText);
                }

                ((JTextField) value).setHorizontalAlignment(JTextField.LEFT);
                value.setSize(dimension);
                value.setPreferredSize(dimension);
                value.setMinimumSize(dimension);
                add(value);
            } else if (valueType.equals(PreferenceType.COLOR)) {

                int valueWidth = (int) (labelWidth * .125);
                dimension = new Dimension(valueWidth, valueHeight);

                value = new ColorPanel();

                // If no value was set but we have a default use it
                if (((valueText == null) || valueText.equals("")) && (defaultValue != null)) {

                    String[] RGB = ((String) defaultValue).split(",");
                    int redValue = new Integer(RGB[0]);
                    int greenValue = new Integer(RGB[1]);
                    int blueValue = new Integer(RGB[2]);
                    Color color = new Color(redValue, greenValue, blueValue);
                    value.setBackground(color);
                } else {

                    String[] RGB = valueText.split(",");
                    int redValue = new Integer(RGB[0]);
                    int greenValue = new Integer(RGB[1]);
                    int blueValue = new Integer(RGB[2]);
                    Color color = new Color(redValue, greenValue, blueValue);
                    value.setBackground(color);
                }

                value.setSize(dimension);
                value.setPreferredSize(dimension);
                value.setMinimumSize(dimension);
                add(value);
            }
            setSize(size);
        }

        /**
         * Method description
         *
         * @return
         */
        public String getExternalPreferenceKey() {
            return externalPreferenceKey;
        }

        /**
         * Method description
         *
         * @return
         */
        public String getKey() {
            return label.getText();
        }

        /**
         * Method description
         *
         * @return
         */
        public String getValue() {

            if (value instanceof JTextField) {
                return ((JTextField) value).getText();
            } else if (value instanceof JCheckBox) {
                boolean isSelected = ((JCheckBox) value).isSelected();
                return isSelected ? "true" : "false";
            } else if (value instanceof ColorPanel) {

                Color selectedColor = ((ColorPanel) value).getBackground();
                if (selectedColor != null) {

                    scratchBuffer.delete(0, scratchBuffer.length());    // Clear
                    int red = selectedColor.getRed();
                    int green = selectedColor.getGreen();
                    int blue = selectedColor.getBlue();
                    scratchBuffer.append(red);
                    scratchBuffer.append(",");
                    scratchBuffer.append(green);
                    scratchBuffer.append(",");
                    scratchBuffer.append(blue);
                }
                return scratchBuffer.toString();
            } else {
                return null;
            }
        }

        /**
         * Method description
         *
         * @param component
         * @return
         */
        @Override
        public Component add(Component component) {

            if (component instanceof JLabel) {
                super.add(component);
            } else if (component instanceof JTextField) {
                super.add(component);
            } else {
                super.add(component);
            }
            return component;
        }


        /*
         * ######################################################################
         * Hiding the super's methods below.
         * ####################################################################
         */

        /**
         * Method description
         *
         * @param component
         * @param arg
         * @return
         */
        @Override
        public Component add(Component component, int arg) {
            return component;
        }

        /**
         * Method description
         *
         * @param arg0
         * @param component
         * @return
         */
        @Override
        public Component add(String arg0, Component component) {
            return component;
        }

        /**
         * Method description
         *
         * @param component
         * @param arg
         */
        @Override
        public void add(Component component, Object arg) {
        }

        /**
         * Method description
         *
         * @param component
         * @param arg
         * @param arg2
         */
        @Override
        public void add(Component component, Object arg, int arg2) {
        }
    }
}
