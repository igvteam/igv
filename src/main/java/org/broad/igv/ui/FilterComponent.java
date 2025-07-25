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

package org.broad.igv.ui;

import org.broad.igv.util.FilterElement;
import org.broad.igv.util.FilterElement.Operator;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class FilterComponent extends javax.swing.JPanel {


    private TrackFilterDialog filterDialog;

    private JComboBox comparisonOperatorComboBox;
    private JComboBox itemComboBox;
    private JButton removeButton;
    private JTextField valueTextField;


    public FilterComponent(TrackFilterDialog filterDialog, List<String> items, FilterElement element) {

        initComponents();

        // Load the Item ComboBox
        this.filterDialog = filterDialog;
        itemComboBox.setModel(new DefaultComboBoxModel(items.toArray()));

        // Load available operators into combobox
        List<String> textForOperators = new ArrayList<String>();
        Operator[] operators = FilterElement.Operator.values();
        for (int i = 0; i < operators.length; i++) {

            // List the operators to skip
            if (operators[i].equals(Operator.GREATER_THAN_OR_EQUAL) ||
                    operators[i].equals(Operator.LESS_THAN_OR_EQUAL)) {
                continue;
            }

            textForOperators.add(operators[i].getValue());
        }
        Collections.sort(textForOperators);
        comparisonOperatorComboBox.setModel(
                new javax.swing.DefaultComboBoxModel(textForOperators.toArray()));


        // Initialize the UI with the FilterElement if present
        if (element != null) {
            itemComboBox.setSelectedItem(element.getAttributeKey());
            valueTextField.setText(element.getValue());
        }

    }


    /**
     * Helper method to convert the string representation of an operator
     * to the appropriate object representation.
     */
    private FilterElement.Operator getOperatorForText(String operatorText) {

        FilterElement.Operator selected = null;

        FilterElement.Operator[] operators = FilterElement.Operator.values();
        for (FilterElement.Operator operator : operators) {

            if (operatorText.equals(operator.getValue())) {
                selected = operator;
                break;
            }
        }

        return selected;
    }

    /**
     * Save the UI content into a non-UI version of the FilterElement
     */
    public FilterElement getFilterElement(boolean matchAll) {
        String attributeKey = itemComboBox.getSelectedItem().toString();
        Operator operator = getOperatorForText((String) comparisonOperatorComboBox.getSelectedItem());
        String expectedValue = valueTextField.getText();
        return new FilterElement(attributeKey, operator, expectedValue);
    }

    public boolean isComplete() {
        return valueTextField.getText() != null && !valueTextField.getText().isEmpty();
    }

    protected void remove() {
        filterDialog.removeComponent(this);
    }


    private void initComponents() {


        setRequestFocusEnabled(false);
        setLayout(new FlowLayout(java.awt.FlowLayout.LEFT, 2, 0));

        itemComboBox = new JComboBox();
        itemComboBox.setMinimumSize(new java.awt.Dimension(50, 27));
        itemComboBox.setPreferredSize(new java.awt.Dimension(150, 27));

        comparisonOperatorComboBox = new JComboBox();
        comparisonOperatorComboBox.setActionCommand("comparisonOperatorComboBoxChanged");
        comparisonOperatorComboBox.setMinimumSize(new java.awt.Dimension(50, 27));
        comparisonOperatorComboBox.setPreferredSize(new java.awt.Dimension(150, 27));

        valueTextField = new JTextField();
        valueTextField.setMaximumSize(new Dimension(32767, 20));
        valueTextField.setMinimumSize(new Dimension(50, 27));
        valueTextField.setPreferredSize(new Dimension(150, 27));

        add(itemComboBox);
        add(comparisonOperatorComboBox);
        add(valueTextField);

        removeButton = new JButton("-");
        removeButton.setFont(new Font("Arial", 0, 14));
        removeButton.setContentAreaFilled(false);
        removeButton.setHorizontalTextPosition(SwingConstants.CENTER);
        removeButton.setPreferredSize(new Dimension(45, 27));
        removeButton.addActionListener(evt -> remove());
        add(removeButton);

    }

}