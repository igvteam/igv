/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

package org.broad.igv.ui;

import com.jidesoft.swing.JideBoxLayout;
import org.apache.log4j.Logger;
import org.broad.igv.track.TrackManager;
import org.broad.igv.ui.panel.*;
import org.broad.igv.ui.util.*;

import javax.swing.*;
import javax.swing.plaf.basic.BasicBorders;
import java.awt.*;
import java.util.*;

import static org.broad.igv.ui.util.UIUtilities.getFileChooser;

/**
 * @author jrobinso
 *
 * Notes;
 *
 * The painting architecture of Swing requires an opaque JComponent to exist in the containment hieararchy above all
 * other components. This is typically provided by way of the content pane. If you replace the content pane, it is
 * recommended that you make the content pane opaque by way of setOpaque(true). Additionally, if the content pane
 * overrides paintComponent, it will need to completely fill in the background in an opaque color in paintComponent.
 *
 * @date Apr 4, 2011
 */
public class IGVContentPane extends JPanel {


    private static Logger log = Logger.getLogger(IGVContentPane.class);

    private JPanel commandBarPanel;
    private IGVCommandBar igvCommandBar;
    private MainPanel mainPanel;
    private ApplicationStatusBar statusBar;

    private TrackManager trackManager;

    /**
     * Creates new form IGV
     */
    public IGVContentPane(TrackManager trackManager) {

        this.trackManager = trackManager;

        // Create components

        setLayout(new BorderLayout());

        commandBarPanel = new JPanel();
        BoxLayout layout = new BoxLayout(commandBarPanel, BoxLayout.PAGE_AXIS);
         
        commandBarPanel.setLayout(layout);
        add(commandBarPanel, BorderLayout.NORTH);

        igvCommandBar = new IGVCommandBar();
        igvCommandBar.setMinimumSize(new Dimension(250, 33));
        igvCommandBar.setBorder(new BasicBorders.MenuBarBorder(Color.GRAY, Color.GRAY));
        igvCommandBar.setAlignmentX(Component.BOTTOM_ALIGNMENT);
        commandBarPanel.add(igvCommandBar);


        mainPanel = new MainPanel(trackManager);
        add(mainPanel, BorderLayout.CENTER);

        statusBar = new ApplicationStatusBar();
        statusBar.setDebugGraphicsOptions(javax.swing.DebugGraphics.NONE_OPTION);
        add(statusBar, BorderLayout.SOUTH);


    }

    public void addCommandBar(JComponent component) {
        component.setBorder(new BasicBorders.MenuBarBorder(Color.GRAY, Color.GRAY));
        component.setAlignmentX(Component.BOTTOM_ALIGNMENT);
        commandBarPanel.add(component);
        commandBarPanel.invalidate();
    }

    public void removeCommandBar(JComponent component)  {
        commandBarPanel.remove(component);
        commandBarPanel.invalidate();
    }

    @Override
    public Dimension getPreferredSize() {
        return UIConstants.preferredSize;
    }

    /**
     * Repaint panels containing data, specifically the dataTrackPanel,
     * featureTrackPanel, and headerPanel.
     */
    public void repaintDataAndHeaderPanels() {
        repaintDataAndHeaderPanels(true);
    }

    public void repaintDataAndHeaderPanels(boolean updateCommandBar) {
        mainPanel.repaint();
        if (updateCommandBar) {
            igvCommandBar.updateCurrentCoordinates();
        }
    }

    public void repaintDataPanels() {
        for (TrackPanelScrollPane tsv : trackManager.getTrackPanelScrollPanes()) {
            tsv.getDataPanel().repaint();
        }

    }

    public void repaintNamePanels() {
        for (TrackPanelScrollPane tsv : trackManager.getTrackPanelScrollPanes()) {
            tsv.getNamePanel().repaint();
        }

    }

    public void repaintStatusAndZoomSlider() {
        igvCommandBar.repaint();
    }



    public void resetFrames() {
        mainPanel.headerPanelContainer.createHeaderPanels();
        for (TrackPanelScrollPane tp : trackManager.getTrackPanelScrollPanes()) {
            tp.getTrackPanel().createDataPanels();
        }

        igvCommandBar.setGeneListMode(FrameManager.isGeneListMode());
        mainPanel.revalidate();
        mainPanel.applicationHeaderPanel.revalidate();
        mainPanel.repaint();
    }



    final public void doRefresh() {

        mainPanel.revalidate();
        repaint();
        //getContentPane().repaint();
    }

    final public void refreshCommandBar() {
        igvCommandBar.updateCurrentCoordinates();
    }


    /**
     * Set the status bar message.  If the message equals "Done." intercept
     * and reset to the default "quite" message,  currently the number of tracks
     * loaded.
     *
     * @param message
     */
    public void setStatusBarMessage(String message) {
        if (message.equals("Done.")) {
            resetStatusMessage();
        }

        ((ApplicationStatusBar) statusBar).setMessage(message);
    }


    /**
     * Add a new data panel set
     */
    public TrackPanelScrollPane addDataPanel(String name) {
        return mainPanel.addDataPanel(name);
    }


    public TrackPanel getDataPanel(String name) {
        TrackPanelScrollPane sp = trackManager.getScrollPane(name);
        if (sp == null) {
            sp = addDataPanel(name);
            trackManager.putScrollPane(name, sp);
        }
        return sp.getTrackPanel();
    }


    public boolean scrollToTrack(String trackName) {
        for (TrackPanelScrollPane sp : trackManager.getTrackPanelScrollPanes()) {
            if (sp.getNamePanel().scrollTo(trackName)) {
                return true;
            }

        }
        return false;
    }

    /**
     * Return an ordered list of track panels.  This method is provided primarily for storing sessions, where
     * the track panels need to be stored in order.
     */
    public java.util.List<TrackPanel> getTrackPanels() {
        return mainPanel.getTrackPanels();
    }



    /**
     * Reset the default status message, which is the number of tracks loaded.
     */
    public void resetStatusMessage() {
        ((ApplicationStatusBar) statusBar).setMessage("" +
                IGV.getInstance().getTrackManager().getVisibleTrackCount() + " tracks loaded");

    }

    public MainPanel getMainPanel() {
        return mainPanel;
    }

    public IGVCommandBar getCommandBar() {
        return igvCommandBar;
    }

    public void chromosomeChanged(String chrName) {
        igvCommandBar.chromosomeChanged(chrName);
    }

    public void updateCurrentCoordinates() {
        igvCommandBar.updateCurrentCoordinates();
    }

    public ApplicationStatusBar getStatusBar() {

        return statusBar;
    }
}
