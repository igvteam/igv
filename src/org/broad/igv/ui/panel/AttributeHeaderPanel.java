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


/*
* TrackPanel.java
*
* Created on Sep 5, 2007, 4:09:39 PM
*
* To change this template, choose Tools | Templates
* and open the template in the editor.
*/

package org.broad.igv.ui.panel;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.PreferenceManager;
import org.broad.igv.track.AttributeManager;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGVMainFrame;

import static org.broad.igv.ui.IGVMainFrame.getInstance;

import org.broad.igv.ui.util.Packable;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import javax.swing.event.MouseInputAdapter;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.geom.AffineTransform;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author jrobinso
 */
public class AttributeHeaderPanel extends JPanel implements Paintable {

    private int attributeCount = 0;
    final static int MAXIMUM_FONT_SIZE = 10;
    final static int ATTRIBUTE_COLUMN_WIDTH = 10;
    final static int COLUMN_BORDER_WIDTH = 2;

    /**
     * Constructs ...
     */
    public AttributeHeaderPanel() {
        setBorder(javax.swing.BorderFactory.createLineBorder(Color.black));
    }


    @Override
    protected void paintComponent(final Graphics graphics) {

        super.paintComponent(graphics);

        // Don't want to destroy the components original graphics
        // context because of the border so we need to create one
        final Graphics2D graphics2 = (Graphics2D) graphics.create();

        final int width = getWidth();
        final int height = getHeight();

        final List<String> keys = getVisibleAttributeKeys();


        final int count = paint(graphics2, keys, width, height);


        for (String key : keys) {

            if (attributeCount != count) {
                addMousableRegion(key);
            }
        }

        //if (attributeCount != count) {
        //    attributeCount = count;
        //    doLayout();
        //}
    }


    public void paintOffscreen(Graphics2D g, Rectangle rect) {

        final int width = rect.width;
        final int height = rect.height;

        final List<String> keys = getVisibleAttributeKeys();

        paint(g, keys, width, height);

    }

    private List<String> getVisibleAttributeKeys() {
        final List<String> keys = AttributeManager.getInstance().getAttributeKeys();
        Set<String> hiddenKeys = AttributeManager.getInstance().getHiddenAttributes();
        keys.removeAll(hiddenKeys);
        return keys;
    }


    private int paint(Graphics2D graphics2, List<String> keys, int width, int height) {
        graphics2.setColor(getBackground());
        graphics2.fillRect(0, 0, width, height);
        graphics2.setColor(Color.BLACK);

        // Divide the remaining space to get column widths
        int columnWidth = getAttributeColumnWidth();

        // Create font and font size
        int fontSize = (int) (0.9 * columnWidth);
        if (fontSize > MAXIMUM_FONT_SIZE) {
            fontSize = MAXIMUM_FONT_SIZE;
        }
        Font font = FontManager.getScalableFont(fontSize);

        final int count = keys.size();
        if (count > 0) {

            if (attributeCount != count) {

                removeAll();    // Remove Mouse Regions
                if (keys != null) {
                    setLayout(new GridLayout(1, count));
                }
            }

            // Change the origin for the text
            AffineTransform transform = AffineTransform.getTranslateInstance(0, height - COLUMN_BORDER_WIDTH);
            graphics2.transform(transform);

            // Now rotate text counter-clockwise 90 degrees
            transform = AffineTransform.getRotateInstance(-Math.PI / 2);
            graphics2.transform(transform);
            graphics2.setFont(font);
            FontMetrics fm = graphics2.getFontMetrics();
            int fontAscent = fm.getHeight();

            int i = 1;
            int x = 0;
            for (String key : keys) {
                int columnLeftEdge = ((COLUMN_BORDER_WIDTH + ATTRIBUTE_COLUMN_WIDTH) * i++);
                x = columnLeftEdge  + ((COLUMN_BORDER_WIDTH + ATTRIBUTE_COLUMN_WIDTH) - fontAscent) / 2;
                graphics2.drawString(key, 0, x);
            }
        }
        return count;
    }


    /**
     * Method description
     *
     * @return
     */
    public int getAttributeColumnWidth() {
        return ATTRIBUTE_COLUMN_WIDTH;
    }

    private static class Region extends JPanel {

        /**
         * Method description
         *
         * @return
         */
        @Override
        public boolean isOpaque() {
            return false;
        }
    }


    private boolean isSortAscending = true;

    private void addMousableRegion(final String tooltip) {

        // Create a clickable region
        final Region region = new Region();
        region.setName(tooltip);
        region.setToolTipText(tooltip);

        MouseInputAdapter listener = new MouseInputAdapter() {

            @Override
            public void mouseClicked(MouseEvent e) {

                String[] selectedAttribute = {region.getName()};
                boolean[] sortAscending = {isSortAscending};
                sortTrackByAttribute(selectedAttribute, sortAscending);
                isSortAscending = !isSortAscending;
            }

            @Override
            public void mouseEntered(MouseEvent e) {

                region.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));

            }

            @Override
            public void mouseExited(MouseEvent e) {
            }

            @Override
            public void mouseReleased(MouseEvent e) {
            }
        };
        region.addMouseMotionListener(listener);
        region.addMouseListener(listener);
        add(region);
    }


    /**
     * Method description
     *
     * @param selectedSortKeys
     * @param isSortAscending
     */
    final public void sortTrackByAttribute(String[] selectedSortKeys, boolean[] isSortAscending) {

        if (selectedSortKeys != null) {

            IGVMainFrame.getInstance().getTrackManager().sortAllTracksByAttributes(selectedSortKeys, isSortAscending);
            IGVMainFrame.getInstance().getContentPane().repaint();
        }
    }

}
