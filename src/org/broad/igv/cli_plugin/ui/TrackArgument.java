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
 * Created by JFormDesigner on Mon Aug 06 15:33:57 EDT 2012
 */

package org.broad.igv.cli_plugin.ui;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import org.broad.igv.cli_plugin.Argument;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.track.DataTrack;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.variant.VariantTrack;

import javax.swing.*;
import java.awt.*;
import java.util.List;

/**
 * @author Jacob Silterra
 */
public class TrackArgument extends ArgumentPanel {
    public TrackArgument(Argument argument) {
        initComponents();
        super.initCommon(argument);

        if (argument != null) {
            List<Track> trackList = IGV.getInstance().getAllTracks();

            Class clazz = getTrackClass(argument);
            Iterable<Track> tracks = Iterables.filter(trackList, clazz);
            List<Track> filteredTrackList = Lists.newArrayList(tracks);
            if(filteredTrackList.size() == 0){
                throw new IllegalStateException("No tracks found of appropriate type; make sure data is loaded");
            }
            trackComboBox.setModel(new DefaultComboBoxModel(filteredTrackList.toArray()));
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
    public Track getValue() {
        return (Track) trackComboBox.getSelectedItem();
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JComboBox trackComboBox;
    // JFormDesigner - End of variables declaration  //GEN-END:variables

    public static class TrackComboBoxRenderer extends DefaultListCellRenderer {

        @Override
        public Component getListCellRendererComponent(JList list, Object value, int index, boolean isSelected, boolean cellHasFocus) {
            Track track = (Track) value;
            String toShow = "No Tracks Found";
            if(track != null){
                toShow = track.getName();
            }
            return super.getListCellRendererComponent(list, toShow, index, isSelected, cellHasFocus);
        }
    }

    public Class getTrackClass(Argument argument) {
        switch (argument.getType()) {
            case VARIANT_TRACK:
                return VariantTrack.class;
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
