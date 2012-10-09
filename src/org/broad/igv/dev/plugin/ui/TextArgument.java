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
 * Created by JFormDesigner on Mon Aug 06 15:14:26 EDT 2012
 */

package org.broad.igv.dev.plugin.ui;

import org.broad.igv.dev.plugin.Argument;

import javax.swing.*;
import java.awt.*;

/**
 * @author User #2
 */
public class TextArgument extends ArgumentPanel {

    public TextArgument() {
        this(null);
    }

    public TextArgument(Argument argument) {
        initComponents();
        super.initCommon(argument);

        if (argument != null) {
            argValue.setText(argument.getDefaultValue());
        }
    }

    @Override
    public String getValue() {
        return argValue.getText();
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        argValue = new JTextField();

        //======== this ========
        setLayout(new BoxLayout(this, BoxLayout.X_AXIS));

        //---- argValue ----
        argValue.setMaximumSize(new Dimension(5000, 28));
        add(argValue);
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JTextField argValue;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
