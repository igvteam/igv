/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */


/*
* To change this template, choose Tools | Templates
* and open the template in the editor.
*/
package org.broad.igv.ui.action;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.commons.math.stat.StatUtils;
import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.List;

/**
 * @author jrobinso
 */
public class SetTrackHeightMenuAction extends MenuAction {

    IGV mainFrame;
    static Logger log = Logger.getLogger(SetTrackHeightMenuAction.class);

    static int lastTrackHeight = -1;

    /**
     * Constructs ...
     *
     * @param label
     * @param mnemonic
     * @param mainFrame
     */
    public SetTrackHeightMenuAction(String label, int mnemonic, IGV mainFrame) {
        super(label, null, mnemonic);
        this.mainFrame = mainFrame;
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

        boolean doRefresh = false;
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

            int status = JOptionPane.showConfirmDialog(mainFrame.getMainFrame(), container, "Set Track Height",
                    JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE, null);

            if ((status == JOptionPane.CANCEL_OPTION) || (status == JOptionPane.CLOSED_OPTION)) {
                return;
            }

            try {
                int newTrackHeight = Integer.parseInt(trackHeightField.getText().trim());
                IGV.getInstance().setAllTrackHeights(newTrackHeight);
                lastTrackHeight = newTrackHeight;
                doRefresh = true;
            }
            catch (NumberFormatException numberFormatException) {
                JOptionPane.showMessageDialog(mainFrame.getMainFrame(), "Track height must be an integer number.");
            }

        }
        finally {

            // Refresh view
            if (doRefresh) {

                // Update the state of the current tracks for drawing purposes

                mainFrame.doRefresh();
            }
            mainFrame.resetStatusMessage();
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
            heights[i] = tracks.get(i).getHeight();
        }
        int medianTrackHeight = (int) Math.round(StatUtils.percentile(heights, 50));
        if (medianTrackHeight > 0) {
            return medianTrackHeight;
        }

        return PreferenceManager.getInstance().getAsInt(PreferenceManager.INITIAL_TRACK_HEIGHT);

    }
}
