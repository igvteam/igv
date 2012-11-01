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

package org.broad.igv.dev.plugin.ui;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import org.broad.igv.dev.plugin.Argument;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.track.DataTrack;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;

import javax.swing.*;
import java.awt.*;
import java.util.List;

/**
 * @author User #2
 */
public class TrackArgument extends ArgumentPanel {
    public TrackArgument(Argument argument) {
        initComponents();
        super.initCommon(argument);

        if (argument != null) {
            List<Track> trackList = IGV.getInstance().getAllTracks();

            Class clazz = getTrackClass(argument);
            Iterable<Track> tracks = Iterables.filter(trackList, clazz);
            trackComboBox.setModel(new DefaultComboBoxModel(Lists.newArrayList(tracks).toArray()));
            trackComboBox.setRenderer(new TrackComboBoxRenderer());
        }
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        trackComboBox = new JComboBox();

        //======== this ========
        setLayout(new BoxLayout(this, BoxLayout.X_AXIS));
        add(trackComboBox);
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    @Override
    public FeatureTrack getValue() {
        return (FeatureTrack) trackComboBox.getSelectedItem();
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
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

    public Class getTrackClass(Argument argument) {
        switch (argument.getType()) {
            case FEATURE_TRACK:
                return FeatureTrack.class;
            case ALIGNMENT_TRACK:
                return AlignmentTrack.class;
            case DATA_TRACK:
                return DataTrack.class;
            default:
                throw new IllegalArgumentException("Argument does not specify a track type; specifies " + argument.getType());
        }
    }

}
