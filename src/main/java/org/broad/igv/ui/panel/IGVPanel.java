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

import org.broad.igv.logging.*;

import javax.swing.*;
import java.awt.*;

/**
 * A panel class that lays out its components in a very specific way so that
 * children (i.e. name, attribute, and data panels) from all instances align vertically
 */
public class IGVPanel extends JPanel implements Paintable {

    Logger log = LogManager.getLogger(IGVPanel.class);

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

            int h = getHeight(); //getPreferredSize().height;
            Component[] children = getComponents();

            int nw = mainPanel.getNamePanelWidth();

            Component namePanel = children[0];
            Component attributePanel = children[1];
            Component dataPanel = children[2];

            namePanel.setBounds(mainPanel.getNamePanelX(), 0, nw, h);
            attributePanel.setBounds(mainPanel.getAttributePanelX(), 0, mainPanel.getAttributePanelWidth(), h);  // Attributes
            dataPanel.setBounds(mainPanel.getDataPanelX(), 0, mainPanel.getDataPanelWidth(), h);

            dataPanel.doLayout();
        }
    }

    public void paintOffscreen(Graphics2D g, Rectangle rect, boolean batch) {

        g.setColor(Color.black);

        Component[] children = getComponents();

        Component namePanel = children[0];
        Rectangle nameRect = new Rectangle(namePanel.getBounds());
        if(nameRect.width > 0) {
            Graphics2D nameGraphics = (Graphics2D) g.create();
            nameGraphics.translate(nameRect.x, 0);
            nameRect.x = 0;
            nameGraphics.setClip(nameRect);
            ((Paintable) namePanel).paintOffscreen(nameGraphics, nameRect, batch);
            nameGraphics.dispose();
        }

        Component attributePanel = children[1];
        Rectangle attRect = new Rectangle(attributePanel.getBounds());
        if(attRect.width > 0) {
            Graphics2D attributeGraphics = (Graphics2D) g.create();
            attributeGraphics.translate(attRect.x, 0);
            attRect.x = 0;
            attributeGraphics.setClip(attRect);
            ((Paintable) attributePanel).paintOffscreen(attributeGraphics, attRect, batch);
            attributeGraphics.dispose();
        }

        Component dataPanel = children[2];
        Rectangle dataRect = new Rectangle(dataPanel.getBounds());
        Graphics2D dataGraphics = (Graphics2D) g.create();
        dataGraphics.translate(dataRect.x, 0);
        dataRect.x = 0;
        g.setClip(dataRect);
        ((Paintable) dataPanel).paintOffscreen(dataGraphics, dataRect, batch);
        dataGraphics.dispose();
    }

    @Override
    public int getSnapshotHeight(boolean batch) {
        return getHeight();
    }
}
