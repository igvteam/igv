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

                    log.debug("Layout: " + frame.getName() + "  x=" + frame.pixelX + "  w=" + frame.getWidthInPixels());

                }
            } catch (Exception e) {
                log.error("Error laying out data panel", e);
            }
        }
    }
}



