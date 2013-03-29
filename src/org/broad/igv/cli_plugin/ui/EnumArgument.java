/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
