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
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.panel;

import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;

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

            log.trace("Layout: " + toString());

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
