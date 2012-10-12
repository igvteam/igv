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
 * Created by JFormDesigner on Wed Aug 08 14:55:42 EDT 2012
 */

package org.broad.igv.dev.plugin.ui;

import com.jidesoft.swing.CheckBoxList;
import org.broad.igv.dev.plugin.Argument;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.TrackWrapper;

import javax.swing.*;
import java.util.ArrayList;
import java.util.List;

/**
 * @author User #2
 */
public class MultiTrackArgument extends ArgumentPanel {

    public MultiTrackArgument(Argument argument) {
        initComponents();
        super.initCommon(argument);

        if (argument != null) {
            trackCheckBoxList.setListData(
                    TrackWrapper.wrapTracks(IGV.getInstance().getFeatureTracks()).toArray()
            );
        }
    }

    @Override
    public List<FeatureTrack> getValue() {
        Object[] rawRet = trackCheckBoxList.getCheckBoxListSelectedValues();
        List<FeatureTrack> trackList = new ArrayList<FeatureTrack>(rawRet.length);
        for (Object obj : rawRet) {
            FeatureTrack track = ((TrackWrapper) obj).getTrack();
            trackList.add(track);
        }
        return trackList;
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        scrollPane1 = new JScrollPane();
        trackCheckBoxList = new CheckBoxList();

        //======== this ========
        setLayout(new BoxLayout(this, BoxLayout.X_AXIS));

        //======== scrollPane1 ========
        {
            scrollPane1.setViewportView(trackCheckBoxList);
        }
        add(scrollPane1);
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JScrollPane scrollPane1;
    private CheckBoxList trackCheckBoxList;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
