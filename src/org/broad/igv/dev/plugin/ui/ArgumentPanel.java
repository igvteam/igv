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
 * Created by JFormDesigner on Mon Aug 06 15:31:39 EDT 2012
 */

package org.broad.igv.dev.plugin.ui;

import org.broad.igv.dev.plugin.Argument;

import javax.swing.*;
import java.awt.*;

/**
 * Base class for an argument row.
 * User: jacob
 * Date: 2012-Aug-08
 */
public class ArgumentPanel extends JPanel {

    public ArgumentPanel() {
        initComponents();
    }

    /**
     * Get the value stored for this argument. Subclasses must override
     *
     * @return
     */
    public Object getValue() {
        return null;
    }

    /**
     * Create the appropriate ArgumentPanel for this argument.
     *
     * @param argument
     * @return
     * @throws IllegalArgumentException If type not found
     */
    public static ArgumentPanel create(Argument argument) {
        ArgumentPanel panel = null;
        switch (argument.getType()) {
            case LONGTEXT:
            case TEXT:
                panel = new TextArgument(argument);
                break;
            case FEATURE_TRACK:
                panel = new TrackArgument(argument);
                break;
            case MULTI_FEATURE_TRACK:
                panel = new MultiTrackArgument(argument);
                break;
            default:
                throw new IllegalArgumentException("Could not create ArgumentPanel for argument of type" + argument.getType());
        }
        return panel;

    }

    protected final void initCommon(Argument argument) {
        if (argument != null) {
            argName.setText(argument.getName() + ":");
        }
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        argName = new JLabel();

        //======== this ========
        setLayout(new BoxLayout(this, BoxLayout.X_AXIS));

        //---- argName ----
        argName.setText("Argument: ");
        argName.setRequestFocusEnabled(false);
        argName.setMinimumSize(new Dimension(80, 16));
        add(argName);
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JLabel argName;
    // JFormDesigner - End of variables declaration  //GEN-END:variables

}
