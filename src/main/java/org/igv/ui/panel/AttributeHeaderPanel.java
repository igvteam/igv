/*
 * TrackPanel.java
 *
 * Created on Sep 5, 2007, 4:09:39 PM
 *
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.igv.ui.panel;

//~--- non-JDK imports --------------------------------------------------------

import org.igv.Globals;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.track.AttributeManager;
import org.igv.ui.FontManager;
import org.igv.ui.IGV;
import org.igv.ui.util.IGVMouseInputAdapter;

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
    private final boolean darkMode;

    Map<String, Boolean> sortOrder = new HashMap();


    public AttributeHeaderPanel() {
        this.darkMode = Globals.isDarkMode();
        setBackground(darkMode ? UIManager.getColor("Panel.background") : new java.awt.Color(255, 255, 255));
        setBorder(javax.swing.BorderFactory.createLineBorder(Color.lightGray));
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

        if(darkMode){
            setBackground(UIManager.getColor("Panel.background"));
        }

        List<String> keys = AttributeManager.getInstance().getVisibleAttributes();

        if (keys != null && keys.size() > 0) {
            final Graphics2D graphics2 = (Graphics2D) graphics.create();
            try {
                if (PreferencesManager.getPreferences().getAntiAliasing()) {
                    graphics2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
                }
                if (darkMode) {
                    graphics2.setColor(Color.white);
                }

                // Divide the remaining space to get column widths
                int columnWidth = getAttributeColumnWidth();

                // Create font and font size
                int fontSize = (int) (0.9 * columnWidth);
                if (fontSize > MAXIMUM_FONT_SIZE) {
                    fontSize = MAXIMUM_FONT_SIZE;
                }
                Font font = FontManager.getFont(fontSize);

                FontMetrics fm = graphics2.getFontMetrics();
                int fontAscent = fm.getHeight();

                // Change the origin for the text
                AffineTransform transform = AffineTransform.getTranslateInstance(0, getHeight() - COLUMN_BORDER_WIDTH);
                graphics2.transform(transform);

                // Now rotate text counter-clockwise 90 degrees
                transform = AffineTransform.getQuadrantRotateInstance(-1);
                graphics2.transform(transform);
                graphics2.setFont(font);

                int i = 1;
                int x;
                for (String key : keys) {
                    int columnLeftEdge = ((COLUMN_BORDER_WIDTH + ATTRIBUTE_COLUMN_WIDTH) * i++);
                    x = columnLeftEdge + ((COLUMN_BORDER_WIDTH + ATTRIBUTE_COLUMN_WIDTH) - fontAscent) / 2;
                    int stringOffset = 2;
                    graphics2.drawString(key, stringOffset, x);
                }
            } finally {
                graphics2.dispose();
            }
        }
    }

    private void addMouseListener() {


        setToolTipText("Click attribute heading to sort");

        MouseInputAdapter listener = new IGVMouseInputAdapter() {

            @Override
            public void igvMouseClicked(MouseEvent e) {

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
            IGV.getInstance().getMainFrame().repaint();
        }
    }

    public void paintOffscreen(Graphics2D g, Rectangle rect, boolean batch) {

        Graphics2D borderGraphics = (Graphics2D) g.create();
        paintComponent(g);
        borderGraphics.setColor(Color.lightGray);
        borderGraphics.drawRect(rect.x, rect.y, rect.width - 1, rect.height - 1);
        borderGraphics.dispose();
    }

    @Override
    public int getSnapshotHeight(boolean batch) {
        return this.getHeight();
    }

}
