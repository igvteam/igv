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

import com.jidesoft.swing.JideScrollPane;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;

/**
 * @author jrobinso
 */
public class TrackPanelScrollPane extends JideScrollPane implements Paintable {

    private static Logger log = Logger.getLogger(TrackPanelScrollPane.class);

    TrackPanel trackPanel;
    boolean isScrolling = false;

    public TrackPanelScrollPane() {
        setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(102, 102, 102)));
        setForeground(new java.awt.Color(153, 153, 153));
        setHorizontalScrollBarPolicy(javax.swing.ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
        setVerticalScrollBarPolicy(javax.swing.ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
        getVerticalScrollBar().setUnitIncrement(16);

        addMouseWheelListener(new MouseWheelListener() {
            public void mouseWheelMoved(MouseWheelEvent mouseWheelEvent) {
                trackPanel.getNamePanel().repaint();
            }
        });

        // A fix for name panel painting problems.   Not sure why this is neccessary, but it is.
        // The adustment listener forces a repaint after a scroll action is complete.
        final JScrollBar sb = this.getVerticalScrollBar();
        sb.addAdjustmentListener(new AdjustmentListener() {
            public void adjustmentValueChanged(AdjustmentEvent adjustmentEvent) {
                if (isScrolling) {
                    if (!adjustmentEvent.getValueIsAdjusting()) {
                        isScrolling = false;
                        trackPanel.getNamePanel().repaint();
                    }
                } else {
                    isScrolling = adjustmentEvent.getValueIsAdjusting();
                }

            }
        });
    }

    @Override
    public void setBackground(Color color) {
        super.setBackground(color);
        if (trackPanel != null)
            trackPanel.setBackground(color);
    }

    @Override
    public void setViewportView(Component trackSetView) {
        if (!(trackSetView instanceof TrackPanel)) {
            throw new IllegalArgumentException("Class TrackPanelScrollPane can only contain a TrackPanel");
        }
        super.setViewportView(trackSetView);
        this.trackPanel = (TrackPanel) trackSetView;
    }

    public TrackPanel getTrackPanel() {
        return trackPanel;
    }


    public String getTrackPanelName() {
        return trackPanel.getName();
    }

    public void minimizeHeight() {
        int prefHeight = trackPanel.getPreferredPanelHeight();
        if (prefHeight < trackPanel.getViewportHeight()) {
            this.setSize(getWidth(), prefHeight);
        }
    }

    public DataPanelContainer getDataPanel() {
        return trackPanel.getDataPanelContainer();
    }

    public boolean isEmpty() {
        return !trackPanel.hasTracks();
    }

    public TrackNamePanel getNamePanel() {
        return trackPanel.getNamePanel();
    }

    @Override
    public void doLayout() {
        log.trace("Layout");
        super.doLayout();    //To change body of overridden methods use File | Settings | File Templates.
        trackPanel.revalidate();
        //trackPanel.doLayout();
    }

    public void paintOffscreen(Graphics2D g, Rectangle tspRect) {

        trackPanel.paintOffscreen(g, tspRect);
    }
}
