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

import org.apache.log4j.Logger;

import javax.swing.*;
import javax.swing.border.Border;
import java.awt.*;

/**
 * @author jrobinso
 * @date Sep 9, 2010
 */
public class DataPanelLayout implements LayoutManager {

    private static Logger log = Logger.getLogger(DataPanelLayout.class);

    int hgap = 5;   // TODO <= make this a function of # of panels

    static Border panelBorder = BorderFactory.createLineBorder(Color.gray);

    public void addLayoutComponent(String s, Component component) {
        // Not used
    }

    public void removeLayoutComponent(Component component) {
        // Not used
    }

    /**
     * Returns the preferred dimensions for this layout given the
     * <i>visible</i> components in the specified target container.
     *
     * @param target the container that needs to be laid out
     *               * @return the preferred dimensions to lay out the
     *               subcomponents of the specified container
     * @see Container
     * @see #minimumLayoutSize
     * @see java.awt.Container#getPreferredSize
     */
    public Dimension preferredLayoutSize(Container target) {
        synchronized (target.getTreeLock()) {
            Dimension dim = new Dimension(0, 0);
            //return dim;
            int nmembers = target.getComponentCount();
            boolean firstVisibleComponent = true;
            for (int i = 0; i < nmembers; i++) {
                Component m = target.getComponent(i);
                if (m.isVisible()) {
                    Dimension d = m.getPreferredSize();
                    dim.height = Math.max(dim.height, d.height);
                    if (firstVisibleComponent) {
                        firstVisibleComponent = false;
                    } else {
                        dim.width += hgap;
                    }
                    dim.width += d.width;
                }
            }
            Insets insets = target.getInsets();
            dim.width += insets.left + insets.right + hgap * 2;
            dim.height += insets.top + insets.bottom;
            return dim;

        }
    }

    public Dimension minimumLayoutSize(Container target) {
        synchronized (target.getTreeLock()) {
            return new Dimension(0, 0);
        }
    }

    public void layoutContainer(Container container) {
        synchronized (container.getTreeLock()) {
            Component[] children = container.getComponents();
            java.util.List<ReferenceFrame> frames = FrameManager.getFrames();
            int h = container.getHeight();

            try {
                //Sometimes the number of children is not the same as the number of frames.
                //Not entirely sure why - Jacob S
                for (int i = 0; i < Math.min(children.length, frames.size()); i++) {

                    Component c = children[i];
                    ReferenceFrame frame = frames.get(i);
                    c.setBounds(frame.pixelX, 0, frame.getWidthInPixels(), h);

                    if (c instanceof JComponent) {
                        if (frame.getWidthInPixels() > 5) {
                            ((JComponent) c).setBorder(panelBorder);
                        } else {
                            ((JComponent) c).setBorder(null);
                        }
                    }

                    log.trace("Layout: " + frame.getName() + "  x=" + frame.pixelX + "  w=" + frame.getWidthInPixels());

                }
            } catch (Exception e) {
                log.error("Error laying out data panel", e);
            }
        }
    }
}



