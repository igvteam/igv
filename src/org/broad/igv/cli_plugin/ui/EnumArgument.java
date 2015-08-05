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
 * Created by JFormDesigner on Fri Mar 29 10:37:25 EDT 2013
 */

package org.broad.igv.cli_plugin.ui;

import org.broad.igv.cli_plugin.Argument;

import javax.swing.*;

/**
 * Argument which can take only one of a fixed set of prescribed values.
 * For instance, a boolean might be yes/no
 * @author jacob
 */
public class EnumArgument extends ArgumentPanel {

    /**
     *
     * @param argument
     * @param values Options the user can choose from
     */
    public EnumArgument(Argument argument, String[] values) {
        initComponents();

        super.initCommon(argument);
        valueComboBox.setModel(new DefaultComboBoxModel(values));

        if (argument != null) {
            valueComboBox.setSelectedItem(argument.getDefaultValue());
        }
    }

    /**
     * @return The string which was selected by the user
     */
    @Override
    public Object getValue() {
        return valueComboBox.getSelectedItem();
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        valueComboBox = new JComboBox();

        //======== this ========
        setLayout(new BoxLayout(this, BoxLayout.X_AXIS));
        add(valueComboBox);
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JComboBox valueComboBox;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
