package org.broad.igv.tools.motiffinder;

import org.broad.igv.feature.CachingFeatureSource;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.MotifTrack;
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
 * Search for a motif (currently can be a regex or IUPAC code)
 */
public class MotifFinderPlugin {

    /**
     * Add menu entry for activating SequenceMatchDialog
     */
    public static JMenuItem getMenuItem() {
        JMenuItem menuItem = new JMenuItem("Find Motif...");
        menuItem.addActionListener(e -> {
            MotifFinderDialog dialog = new MotifFinderDialog(IGV.getInstance().getMainFrame());
            dialog.setVisible(true);
            handleDialogResult(dialog);
        });

        return menuItem;
    }

    static void handleDialogResult(MotifFinderDialog dialog) {
        String[] pattern = dialog.getInputPattern();
        if (pattern != null) {
            String[] posTrackName = dialog.getPosTrackName();
            String[] negTrackName = dialog.getNegTrackName();
            addTracksForPatterns(pattern, posTrackName, negTrackName);
        }
    }

    /**
     * Generate motif-finding track and add it to IGV
     *
     * @param pattern
     * @param posTrackNames
     * @param negTrackNames
     * @return
     */
    static List<Track> addTracksForPatterns(String[] pattern, String[] posTrackNames, String[] negTrackNames) {
        List<Track> trackList = generateTracksForPatterns(pattern, posTrackNames, negTrackNames);
        IGV.getInstance().addTracks(trackList);
        return trackList;
    }

    /**
     * Generate motif-finding tracks for the given pattern, do not add them to anything
     *
     * @param patterns
     * @param posTrackNames
     * @param negTrackNames
     * @return
     */
    static List<Track> generateTracksForPatterns(String[] patterns, String[] posTrackNames, String[] negTrackNames) {

        Color[] colors = {null, Color.RED};
        Strand[] strands = {Strand.POSITIVE, Strand.NEGATIVE};
        List<Track> trackList = new ArrayList<Track>(2 * posTrackNames.length);

        if (patterns != null) {
            for (int pi = 0; pi < patterns.length; pi++) {
                String pattern = patterns[pi];
                String[] curTrackNames = new String[]{posTrackNames[pi], negTrackNames[pi]};
                for (int ci = 0; ci < curTrackNames.length; ci++) {
                    String tName = curTrackNames[ci];
                    if (tName == null) continue;

                    MotifTrack track = new MotifTrack(tName, pattern, strands[ci]);
                    if (colors[ci] != null) track.setColor(colors[ci]);

                    track.setDisplayMode(Track.DisplayMode.SQUISHED);
                    trackList.add(track);
                }
            }

        }
        return trackList;
    }

    public String run(List<String> args) {
        String cmd = args.get(0);
        if (cmd.equalsIgnoreCase("find")) {
            String pattern = args.get(1);
            String[] patterns = new String[]{pattern};
            String shrtPattern = StringUtils.checkLength(pattern, MotifFinderDialog.MaxTrackNameLength);
            String[] posName = new String[]{shrtPattern + " Positive"};
            String[] negName = new String[]{shrtPattern + " Negative"};
            addTracksForPatterns(patterns, posName, negName);
            return "OK";
        } else {
            return "ERROR: Unknown command " + cmd + " for plugin " + getClass().getName();
        }
    }
}
