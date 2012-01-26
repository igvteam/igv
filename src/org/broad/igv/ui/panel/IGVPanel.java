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
package org.broad.igv.ui.panel;

import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.util.Arrays;
import java.util.Comparator;

/**
 * A panel class that lays out its components in a very specific way so that
 * children (i.e. name, attribute, and data panels) from all instances align vertically
 */
public class IGVPanel extends JPanel implements Paintable {

    Logger log = Logger.getLogger(IGVPanel.class);

    // Backpointer to parent
    MainPanel mainPanel;

    public IGVPanel(MainPanel mainPanel) {
        setLayout(null); //new TrackPanelLayout());
        this.mainPanel = mainPanel;
    }

    public int getViewportHeight() {

        Container parent = getParent();
        return parent == null ? 0 : parent.getHeight();
    }

    public TrackPanelScrollPane getScrollPane() {

        TrackPanelScrollPane scollpane = null;
        Container parent = getParent();

        if (parent instanceof JViewport) {
            scollpane = (TrackPanelScrollPane) ((JViewport) parent).getParent();
        }
        return scollpane;
    }

    @Override
    public void doLayout() {
        synchronized (getTreeLock()) {

            log.debug("Layout: " + toString());

            int h = getHeight(); //getPreferredSize().height;
            Component[] children = getComponents();

            int nw = mainPanel.getNamePanelWidth();
            int grabBarWidth = nw - 10;

            int idx = 0;
            if (children.length > 3) {
                children[idx++].setBounds(mainPanel.getNamePanelX(), 0, 10, h);
                children[idx++].setBounds(mainPanel.getNamePanelX() + 10, 0, nw - 10, h);
            } else {
                children[idx++].setBounds(mainPanel.getNamePanelX(), 0, nw, h);

            }
            children[idx++].setBounds(mainPanel.getAttributePanelX(), 0, mainPanel.getAttributePanelWidth(), h);
            children[idx].setBounds(mainPanel.getDataPanelX(), 0, mainPanel.getDataPanelWidth(), h);

            children[idx].doLayout();
        }
    }

    public void paintOffscreen(Graphics2D g, Rectangle rect) {

        int h = rect.height;

        Component[] children = getComponents();
        // name panel starts at offset=0

        g.translate(mainPanel.getNamePanelX(), 0);

        Rectangle nameRect = new Rectangle(children[0].getBounds());
        nameRect.height = h;
        g.setClip(nameRect);
        ((Paintable) children[0]).paintOffscreen(g, nameRect);

        int dx = mainPanel.getAttributePanelX() - mainPanel.getNamePanelX();
        g.translate(dx, 0);
        Rectangle attRect = new Rectangle(0, 0, children[1].getWidth(), h);
        g.setClip(attRect);
        ((Paintable) children[1]).paintOffscreen(g, attRect);

        dx = mainPanel.getDataPanelX() - mainPanel.getAttributePanelX();
        g.translate(dx, 0);
        Rectangle dataRect = new Rectangle(0, 0, mainPanel.getDataPanelWidth(), h);
        g.setClip(dataRect);
        ((Paintable) children[2]).paintOffscreen(g, dataRect);



    }


}
