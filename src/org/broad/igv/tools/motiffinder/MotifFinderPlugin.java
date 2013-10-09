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

package org.broad.igv.tools.motiffinder;

import org.broad.igv.dev.api.IGVPlugin;
import org.broad.igv.dev.api.batch.Command;
import org.broad.igv.feature.CachingFeatureSource;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.PanelName;
import org.broad.igv.util.StringUtils;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.List;

/**
 * Plugin for searching for a motif (currently can be a regex or IUPAC code)
 * This plugin provides for a dialog so that the user can enter motifs manually,
 * also a batch {@link org.broad.igv.dev.api.batch.Command}.
* @author jacob
* @date 2013-Oct-09
*/
public class MotifFinderPlugin implements IGVPlugin, Command {

    /**
     * Add menu entry for activating SequenceMatchDialog
     */
    @Override
    public void init() {
        JMenuItem menuItem = new JMenuItem("Find Motif...");
        menuItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                MotifFinderDialog dialog = new MotifFinderDialog(IGV.getMainFrame());

                dialog.setVisible(true);

                String pattern = dialog.getInputPattern();

                if (pattern != null) {
                    String posTrackName = dialog.getPosTrackName();
                    String negTrackName = dialog.getNegTrackName();
                    addTracksForPattern(pattern, posTrackName, negTrackName);
                }
            }
        });

        IGV.getInstance().addOtherToolMenu(menuItem);
    }

    /**
     * Generate motif-finding track and add it to IGV
     * @param pattern
     * @param posTrackName
     * @param negTrackName
     * @return
     */
    static List<Track> addTracksForPattern(String pattern, String posTrackName, String negTrackName){
        List<Track> trackList = generateTracksForPattern(pattern, posTrackName, negTrackName);
        IGV.getInstance().addTracks(trackList, PanelName.FEATURE_PANEL);
        return trackList;
    }

    /**
     * Generate motif-finding tracks for the given pattern, do not add them to anything
     * @param pattern
     * @param posTrackName
     * @param negTrackName
     * @return
     */
    private static List<Track> generateTracksForPattern(String pattern, String posTrackName, String negTrackName){

        String[] trackNames = {posTrackName, negTrackName};
        Color[] colors = {null, Color.RED};
        Strand[] strands = {Strand.POSITIVE, Strand.NEGATIVE};
        List<Track> trackList = new ArrayList<Track>(trackNames.length);

        if (pattern != null) {
            for(int ii=0; ii < trackNames.length; ii++){
                String tName = trackNames[ii];
                if(tName == null) continue;

                MotifFinderSource src = new MotifFinderSource(pattern, strands[ii], GenomeManager.getInstance().getCurrentGenome());
                CachingFeatureSource cachingSrc= new CachingFeatureSource(src);

                FeatureTrack track = new FeatureTrack(tName, tName, cachingSrc);
                if(colors[ii] != null) track.setColor(colors[ii]);

                track.setDisplayMode(Track.DisplayMode.SQUISHED);
                trackList.add(track);
            }
        }
        return trackList;
    }

    @Override
    public String run(List<String> args) {
        String cmd = args.get(0);
        if(cmd.equalsIgnoreCase("find")){
            String pattern = args.get(1);
            String shrtPattern = StringUtils.checkLength(pattern, MotifFinderDialog.MaxTrackNameLength);
            String posName = shrtPattern + " Positive";
            String negName = shrtPattern + " Negative";
            addTracksForPattern(pattern, posName, negName);
            return "OK";
        }else{
            return "ERROR: Unknown command " + cmd + " for plugin " + getClass().getName();
        }
    }
}
