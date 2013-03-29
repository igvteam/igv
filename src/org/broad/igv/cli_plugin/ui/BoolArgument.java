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
 * Created by JFormDesigner on Fri Mar 29 11:33:04 EDT 2013
 */

package org.broad.igv.cli_plugin.ui;

import org.broad.igv.cli_plugin.Argument;

import javax.swing.*;

/**
 * UI representation of a boolean value; shows a checkbox. Checked = true
 * @author jacob
 */
public class BoolArgument extends ArgumentPanel {

    public BoolArgument(Argument argument) {
        initComponents();
        super.initCommon(argument);

        if(argument != null){
            valueCheckBox.setSelected(Boolean.parseBoolean(argument.getDefaultValue()));
        }
    }

    @Override
    public Object getValue() {
        return valueCheckBox.isSelected();
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        valueCheckBox = new JCheckBox();

        //======== this ========
        setLayout(new BoxLayout(this, BoxLayout.X_AXIS));
        add(valueCheckBox);
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JCheckBox valueCheckBox;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
