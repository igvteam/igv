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
 * Created by JFormDesigner on Mon Aug 06 15:31:39 EDT 2012
 */

package org.broad.igv.cli_plugin.ui;

import org.broad.igv.cli_plugin.Argument;

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

    void setArgName(String newName){
        this.argName.setText(newName);
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
            case LOCUS:
                break;
            case BOOL:
                panel = new BoolArgument(argument);
                break;
            case LONGTEXT:
            case TEXT:
                panel = new TextArgument(argument);
                break;
            case VARIANT_TRACK:
            case DATA_TRACK:
            case ALIGNMENT_TRACK:
            case FEATURE_TRACK:
                panel = new TrackArgument(argument);
                break;
            case MULTI_FEATURE_TRACK:
                panel = new MultiTrackArgument(argument);
                break;
            default:
                throw new IllegalArgumentException("Could not create ArgumentPanel for argument of type " + argument.getType());
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
