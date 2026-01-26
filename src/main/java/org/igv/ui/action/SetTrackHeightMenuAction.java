package org.igv.ui.action;


import org.apache.commons.math3.stat.StatUtils;
import org.igv.logging.*;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.track.Track;
import org.igv.ui.IGV;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.List;

/**
 * @author jrobinso
 */
public class SetTrackHeightMenuAction extends MenuAction {

    IGV igv;
    static Logger log = LogManager.getLogger(SetTrackHeightMenuAction.class);

    static int lastTrackHeight = -1;

    /**
     * Constructs ...
     *
     * @param label
     * @param mnemonic
     * @param igv
     */
    public SetTrackHeightMenuAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
    }

    /**
     * Method description
     *
     * @param e
     */
    @Override
    public void actionPerformed(ActionEvent e) {
        doSetTrackHeight();

    }

    /**
     * Method description
     */
    final public void doSetTrackHeight() {

        boolean repaint = false;
        try {
            JPanel container = new JPanel();
            JLabel trackHeightLabel = new JLabel("Track Height (pixels)");
            JTextField trackHeightField = new JTextField();
            Dimension preferredSize = trackHeightField.getPreferredSize();
            trackHeightField.setPreferredSize(new Dimension(50, (int) preferredSize.getHeight()));
            container.add(trackHeightLabel);
            container.add(trackHeightField);

            int repTrackHeight = getRepresentativeTrackHeight();
            trackHeightField.setText(String.valueOf(repTrackHeight));

            int status = JOptionPane.showConfirmDialog(igv.getMainFrame(), container, "Set Track Height",
                    JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE, null);

            if ((status == JOptionPane.CANCEL_OPTION) || (status == JOptionPane.CLOSED_OPTION)) {
                return;
            }

            try {
                int newTrackHeight = Integer.parseInt(trackHeightField.getText().trim());
                IGV.getInstance().setAllTrackHeights(newTrackHeight);
                lastTrackHeight = newTrackHeight;
                repaint = true;
            }
            catch (NumberFormatException numberFormatException) {
                JOptionPane.showMessageDialog(igv.getMainFrame(), "Track height must be an integer number.");
            }

        }
        finally {

            // Refresh view
            if (repaint) {

                // Update the state of the current tracks for drawing purposes

                igv.repaint();
            }
            igv.resetStatusMessage();
        }

    }

    /**
     * Return a representative track height to use as the default.  For now
     * using the median track height.
     *
     * @return
     */
    private int getRepresentativeTrackHeight() {

        if (lastTrackHeight > 0) {
            return lastTrackHeight;
        }

        // Get all tracks except the gene track
        List<Track> tracks = IGV.getInstance().getAllTracks();


        double[] heights = new double[tracks.size()];
        for (int i = 0; i < tracks.size(); i++) {
            heights[i] = tracks.get(i).getContentHeight();
        }
        int medianTrackHeight = (int) Math.round(StatUtils.percentile(heights, 50));
        if (medianTrackHeight > 0) {
            return medianTrackHeight;
        }

        return PreferencesManager.getPreferences().getAsInt(Constants.INITIAL_TRACK_HEIGHT);

    }
}
