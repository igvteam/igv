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

package org.broad.igv.ui;

import org.apache.log4j.Logger;
import org.broad.igv.ui.panel.MainPanel;
import org.broad.igv.ui.panel.TrackPanel;
import org.broad.igv.ui.util.ApplicationStatusBar;

import javax.swing.*;
import javax.swing.plaf.basic.BasicBorders;
import java.awt.*;


/**
 * @author jrobinso
 *         <p/>
 *         Notes;
 *         <p/>
 *         The painting architecture of Swing requires an opaque JComponent to exist in the containment hieararchy above all
 *         other components. This is typically provided by way of the content pane. If you replace the content pane, it is
 *         recommended that you make the content pane opaque by way of setOpaque(true). Additionally, if the content pane
 *         overrides paintComponent, it will need to completely fill in the background in an opaque color in paintComponent.
 * @date Apr 4, 2011
 */
public class IGVContentPane extends JPanel {


    private static Logger log = Logger.getLogger(IGVContentPane.class);

    private JPanel commandBarPanel;
    private IGVCommandBar igvCommandBar;
    private MainPanel mainPanel;
    private ApplicationStatusBar statusBar;

    private IGV igv;

    /**
     * Creates new form IGV
     */
    public IGVContentPane(IGV igv) {

        this.igv = igv;

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


        mainPanel = new MainPanel(igv);
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

    public void removeCommandBar(JComponent component) {
        commandBarPanel.remove(component);
        commandBarPanel.invalidate();
    }

    @Override
    public Dimension getPreferredSize() {
        return UIConstants.preferredSize;
    }


    public void repaintDataPanels() {
        for (TrackPanel tp : mainPanel.getTrackPanels()) {
            tp.getScrollPane().getDataPanel().repaint();
        }
    }

    public void revalidateDataPanels() {
        for (TrackPanel tp : mainPanel.getTrackPanels()) {
            tp.getScrollPane().getDataPanel().revalidate();
        }
    }


    final public void doRefresh() {

        mainPanel.revalidate();
        repaint();
        //getContentPane().repaint();
    }

    /**
     * Reset the default status message, which is the number of tracks loaded.
     */
    public void resetStatusMessage() {
        statusBar.setMessage("" + igv.getVisibleTrackCount() + " tracks loaded");

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
