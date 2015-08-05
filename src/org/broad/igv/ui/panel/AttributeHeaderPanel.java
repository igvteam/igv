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
* TrackPanel.java
*
* Created on Sep 5, 2007, 4:09:39 PM
*
* To change this template, choose Tools | Templates
* and open the template in the editor.
*/

package org.broad.igv.ui.panel;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.track.AttributeManager;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;

import javax.swing.*;
import javax.swing.event.MouseInputAdapter;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.geom.AffineTransform;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class AttributeHeaderPanel extends JPanel implements Paintable {

    final static int MAXIMUM_FONT_SIZE = 10;
    public final static int ATTRIBUTE_COLUMN_WIDTH = 10;
    public final static int COLUMN_BORDER_WIDTH = 1;

    Map<String, Boolean> sortOrder = new HashMap();


    public AttributeHeaderPanel() {
        setBorder(javax.swing.BorderFactory.createLineBorder(Color.black));
        setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
        addMouseListener();
    }


    public int getAttributeColumnWidth() {
        return ATTRIBUTE_COLUMN_WIDTH;
    }


    private String getAttributeHeading(int x) {
        int idx = x / (ATTRIBUTE_COLUMN_WIDTH + COLUMN_BORDER_WIDTH);
        List<String> keys = AttributeManager.getInstance().getVisibleAttributes();
        if (idx < keys.size()) {
            return keys.get(idx);
        } else {
            return null;
        }
    }

    @Override
    protected void paintComponent(final Graphics graphics) {

        super.paintComponent(graphics);


        List<String> keys = AttributeManager.getInstance().getVisibleAttributes();

        if (keys != null && keys.size() > 0) {

            final Graphics2D graphics2 = (Graphics2D) graphics.create();
            if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.ENABLE_ANTIALISING)) {
                graphics2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
            }

            // Divide the remaining space to get column widths
            int columnWidth = getAttributeColumnWidth();

            // Create font and font size
            int fontSize = (int) (0.9 * columnWidth);
            if (fontSize > MAXIMUM_FONT_SIZE) {
                fontSize = MAXIMUM_FONT_SIZE;
            }
            Font font = FontManager.getFont(fontSize);

            // Change the origin for the text
            AffineTransform transform = AffineTransform.getTranslateInstance(0, getHeight() - COLUMN_BORDER_WIDTH);
            graphics2.transform(transform);

            // Now rotate text counter-clockwise 90 degrees
            transform = AffineTransform.getQuadrantRotateInstance(-1);
            graphics2.transform(transform);
            graphics2.setFont(font);
            graphics2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
            FontMetrics fm = graphics2.getFontMetrics();
            int fontAscent = fm.getHeight();

            int i = 1;
            int x;
            for (String key : keys) {
                int columnLeftEdge = ((COLUMN_BORDER_WIDTH + ATTRIBUTE_COLUMN_WIDTH) * i++);
                x = columnLeftEdge + ((COLUMN_BORDER_WIDTH + ATTRIBUTE_COLUMN_WIDTH) - fontAscent) / 2;
                String toDraw = key;
                int stringOffset = 2;
                graphics2.drawString(toDraw, stringOffset, x);
            }
        }
    }

    private void addMouseListener() {


        setToolTipText("Click attribute heading to sort");

        MouseInputAdapter listener = new MouseInputAdapter() {

            @Override
            public void mouseClicked(MouseEvent e) {

                String attKey = getAttributeHeading(e.getX());
                if (attKey != null) {
                    Boolean tmp = sortOrder.get(attKey);
                    boolean sortAscending = tmp == null ? true : tmp.booleanValue();
                    sortTrackByAttribute(attKey, sortAscending);
                    sortOrder.put(attKey, !sortAscending);
                }
            }

            @Override
            public void mouseMoved(MouseEvent e) {
                String attKey = getAttributeHeading(e.getX());
                if (attKey != null) {
                    AttributeManager.ColumnMetaData md = AttributeManager.getInstance().getColumnMetaData(attKey);
                    if (md != null) {
                        StringBuffer buf = new StringBuffer("<html>" + attKey + "<br>Click to sort");
                        //buf.append("<br>Is numeric = " + md.isNumeric() + "<br>");
                        //buf.append("Is diverging = " + md.isDiverging() + "<br>");
                        //buf.append("getUniqueCount = " + md.getUniqueCount() + "<br>");
                        //buf.append("getTotalCount = " + md.getTotalCount() + "<br>");
                        //buf.append("getUniqueRatio = " + md.getUniqueRatio() + "<br>");
                        // buf.append("# unique = " + md.uniqueValues.size() + "<br>");
                        // buf.append("# total = " + md.totalCount + "<br>");
                        // buf.append("# numeric = " + md.numericCount + "<br>");
                        setToolTipText(buf.toString());
                        return;
                    }
                    setToolTipText("Click attribute heading to sort");
                }

            }

        };
        addMouseMotionListener(listener);
        addMouseListener(listener);
    }


    final public void sortTrackByAttribute(String sortKey, boolean isSortAscending) {

        if (sortKey != null) {
            IGV.getInstance().sortAllTracksByAttributes(new String[]{sortKey}, new boolean[]{isSortAscending});
            IGV.getMainFrame().repaint();
        }
    }

    public void paintOffscreen(Graphics2D g, Rectangle rect) {
        paintComponent(g);

        Color c = g.getColor();
        g.setColor(Color.darkGray);
        g.drawRect(rect.x, rect.y, rect.width, rect.height);
        g.setColor(c);            //super.paintBorder(g);

    }
}
