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

package org.broad.igv.ui.panel;

import org.broad.igv.lists.GeneList;
import org.broad.igv.ui.IGV;

import javax.swing.*;
import javax.swing.border.Border;
import java.awt.*;

/**
 * @author Jim Robinson
 * @date 5/25/12
 */
public class PanelLayoutManager {

//    HeaderPanelContainer headerPanelContainer;
//    DataPanelContainer dataPanelContainer;
//
//    public void layoutContainer(Container container) {
//        synchronized (container.getTreeLock()) {
//            Component[] children = container.getComponents();
//            //java.util.List<ReferenceFrame> frames = FrameManager.getFrames();
//            int h = container.getHeight();
//
//            try {
//                //Sometimes the number of children is not the same as the number of frames.
//                //Not entirely sure why - Jacob S
//                for (int i = 0; i < Math.min(children.length, frames.size()); i++) {
//
//                    Component c = children[i];
//                    ReferenceFrame frame = frames.get(i);
//                    c.setBounds(frame.pixelX, 0, frame.getWidthInPixels(), h);
//
//                    if (c instanceof JComponent) {
//                        if (frame.getWidthInPixels() > 5) {
//                            ((JComponent) c).setBorder(panelBorder);
//                        } else {
//                            ((JComponent) c).setBorder(null);
//                        }
//                    }
//
//                    log.debug("Layout: " + frame.getName() + "  x=" + frame.pixelX + "  w=" + frame.getWidthInPixels());
//
//                }
//            } catch (Exception e) {
//                log.error("Error laying out data panel", e);
//            }
//        }
//    }
//
//    public void createHeaderPanels() {
//
//        removeAll();
//        frames = FrameManager.getFrames();
//
//        contentPanel = new JPanel();
//        contentPanel.setLayout(new DataPanelLayout());
//        for (ReferenceFrame f : frames) {
//            HeaderPanel dp = new HeaderPanel(f);
//            dp.setBackground(getBackground());
//            contentPanel.add(dp);
//        }
//        add(contentPanel, BorderLayout.CENTER);
//
//        if (FrameManager.isGeneListMode()) {
//            GeneList gl = IGV.getInstance().getSession().getCurrentGeneList();
//            String name = gl.getDisplayName();
//            JLabel label = new JLabel(name, JLabel.CENTER);
//            Border border = BorderFactory.createLineBorder(Color.lightGray);
//            label.setBorder(border);
//            add(label, BorderLayout.NORTH);
//        }
//
//        invalidate();
//    }
//
//    public void createDataPanels() {
//        removeAll();
//        for (ReferenceFrame f : FrameManager.getFrames()) {
//            DataPanel dp = new DataPanel(f, this);
//            add(dp);
//        }
//        invalidate();
//    }

}
