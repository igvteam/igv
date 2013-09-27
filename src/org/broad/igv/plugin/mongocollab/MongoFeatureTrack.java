/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.plugin.mongocollab;

import com.mongodb.DBCollection;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackClickEvent;
import org.broad.igv.track.TrackMenuUtils;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.tribble.Feature;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Arrays;

/**
 * @author jacob
 * @date 2013-Sep-27
 */
public class MongoFeatureTrack extends FeatureTrack{

    public MongoFeatureTrack(String id, String name, MongoFeatureSource source) {
        super(id, name, source);
    }

    private JMenuItem createAddAnnotMenuEntry(TrackClickEvent te){
        ReferenceFrame frame = te.getFrame();
        boolean hasFrame = frame != null;
        if(!hasFrame) return null;

        final String chr = frame.getChrName();
        final int start = (int) te.getChromosomePosition();
        final int end = (int) Math.ceil(frame.getChromosomePosition(te.getMouseEvent().getX() + 1));

        JMenuItem item = new JMenuItem("Add annotation to " + MongoFeatureTrack.this.getName());
        item.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                Feature feat = new BasicFeature(chr, start, end);
                FeatureAnnotDialog dialog = new FeatureAnnotDialog(IGV.getMainFrame(), getCollection(), feat);
                dialog.setVisible(true);
            }
        });
        return item;
    }

    private DBCollection getCollection(){
        return ((MongoFeatureSource) source).getCollection();
    }

    @Override
    public IGVPopupMenu getPopupMenu(TrackClickEvent te) {
        String title = getName();
        IGVPopupMenu menu = TrackMenuUtils.getPopupMenu(Arrays.<Track>asList(this), title, te);
        menu.add(createAddAnnotMenuEntry(te));
        return menu;
    }
}
