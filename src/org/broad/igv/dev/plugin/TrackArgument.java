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
 * Created by JFormDesigner on Mon Aug 06 15:33:57 EDT 2012
 */

package org.broad.igv.dev.plugin;

import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;

import javax.swing.*;
import java.awt.*;

/**
 * @author User #2
 */
public class TrackArgument extends ArgumentPanel {
    public TrackArgument(Argument argument) {
        initComponents();

        if (argument != null) {
            argName.setText(argument.getName() + ":");
            cmdArg.setText(argument.getCmdArg());
            trackComboBox.setModel(new DefaultComboBoxModel((IGV.getInstance().getFeatureTracks()).toArray()));
            trackComboBox.setRenderer(new TrackComboBoxRenderer());
        }
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        argName = new JLabel();
        cmdArg = new JLabel();
        trackComboBox = new JComboBox();

        //======== this ========
        setLayout(new BoxLayout(this, BoxLayout.X_AXIS));

        //---- argName ----
        argName.setText("Argument: ");
        add(argName);

        //---- cmdArg ----
        cmdArg.setText("text");
        add(cmdArg);
        add(trackComboBox);
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    @Override
    public FeatureTrack getValue() {
        return (FeatureTrack) trackComboBox.getSelectedItem();
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JLabel argName;
    private JLabel cmdArg;
    private JComboBox trackComboBox;
    // JFormDesigner - End of variables declaration  //GEN-END:variables

    public static class TrackComboBoxRenderer extends DefaultListCellRenderer {

        @Override
        public Component getListCellRendererComponent(JList list, Object value, int index, boolean isSelected, boolean cellHasFocus) {
            Track track = (Track) value;
            String toShow = track.getName();
            return super.getListCellRendererComponent(list, toShow, index, isSelected, cellHasFocus);
        }
    }
}
