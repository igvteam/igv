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
 * Created by JFormDesigner on Wed Aug 08 14:55:42 EDT 2012
 */

package org.broad.igv.cli_plugin.ui;

import com.jidesoft.swing.CheckBoxList;
import org.broad.igv.cli_plugin.Argument;
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
            FeatureTrack track = (FeatureTrack) ((TrackWrapper) obj).getTrack();
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
