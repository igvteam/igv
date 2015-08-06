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

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.track.AttributeManager;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackClickEvent;
import org.broad.igv.track.TrackGroup;
import org.broad.igv.ui.IGV;

import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.util.Packable;

import javax.swing.event.MouseInputAdapter;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class AttributePanel extends TrackPanelComponent implements Packable, Paintable {

    private static Logger log = Logger.getLogger(AttributePanel.class);


    /**
     * Constructs ...
     */
    public AttributePanel(TrackPanel trackPanel) {
        super(trackPanel);
        setBackground(Color.lightGray);
        setBorder(javax.swing.BorderFactory.createLineBorder(Color.black));
        init();
    }


    @Override
    protected void paintComponent(Graphics g) {

        super.paintComponent(g);
        Rectangle visibleRect = getVisibleRect();
        removeMousableRegions();
        paintOffscreen((Graphics2D) g, visibleRect);

    }


    public void paintOffscreen(Graphics2D g, Rectangle rect) {

        List<String> names = AttributeManager.getInstance().getAttributeNames();

        if (names == null) {
            return;
        }

        final Set<String> hiddenAttributes = IGV.getInstance().getSession().getHiddenAttributes();
        if (hiddenAttributes != null) names.removeAll(hiddenAttributes);

        if (names.size() > 0) {


            // Get the current tracks
            TrackPanel trackPanel = (TrackPanel) getParent();
            Collection<TrackGroup> groups = trackPanel.getGroups();


            if (!groups.isEmpty()) {

                // int attributeColumnWidth = getAttributeColumnWidth();
                final Graphics2D graphics2D = (Graphics2D) g.create();
                graphics2D.setColor(Color.BLACK);

                final Graphics2D greyGraphics = (Graphics2D) g.create();
                greyGraphics.setColor(UIConstants.LIGHT_GREY);

                final Graphics2D borderGraphics = (Graphics2D) g.create();
                borderGraphics.setColor(Color.lightGray);

                final int left = AttributeHeaderPanel.COLUMN_BORDER_WIDTH;
                int regionY = 0;
                final int bottom = rect.y + rect.height;

                for (Iterator<TrackGroup> groupIter = groups.iterator(); groupIter.hasNext(); ) {
                    TrackGroup group = groupIter.next();

                    // Out of view?
                    if (regionY > bottom) {
                        break;
                    }

                    if (group.isVisible()) {
                        if (groups.size() > 1) {
                            greyGraphics.fillRect(0, regionY + 1, getWidth(), UIConstants.groupGap - 1);
                            regionY += UIConstants.groupGap;
                        }


                        if (group.isDrawBorder()) {
                            g.drawLine(0, regionY - 1, getWidth(), regionY - 1);
                        }

                        for (Track track : group.getTracks()) {
                            if (track == null) continue;
                            int trackHeight = track.getHeight();
                            if (regionY > bottom) {
                                break;
                            }

                            if (track.isVisible()) {
                                int border = trackHeight < 5 ? 0 : 1;
                                if (regionY + trackHeight >= rect.y) {
                                    Rectangle trackRectangle = new Rectangle(left, regionY+border, getWidth(), trackHeight-border);
                                    track.renderAttributes(graphics2D, trackRectangle, rect, names, mouseRegions);
                                    //regionY = draw(names, track, regionX, regionY, attributeColumnWidth, track.getHeight(), graphics2D);
                                }
                                regionY += trackHeight;
                            }
                        }

                        if (group.isDrawBorder()) {
                            g.drawLine(0, regionY, getWidth(), regionY);
                        }
                    }
                }

                // Border

                final int colWidth = AttributeHeaderPanel.ATTRIBUTE_COLUMN_WIDTH + AttributeHeaderPanel.COLUMN_BORDER_WIDTH;
                for (int x = 1; x < rect.x + rect.width; x += colWidth) {
                    borderGraphics.fillRect(x + AttributeHeaderPanel.ATTRIBUTE_COLUMN_WIDTH, rect.y,
                            AttributeHeaderPanel.COLUMN_BORDER_WIDTH, rect.height);
                }

                graphics2D.dispose();
                greyGraphics.dispose();
                borderGraphics.dispose();

            }

        }

    }


    /*private int draw(List<String> names, Track track, int trackX, int trackY, int trackWidth, int trackHeight,
                     Graphics2D graphics) {

        for (String name : names) {

            String key = name.toUpperCase();
            String attributeValue = track.getAttributeValue(key);

            if (attributeValue != null) {

                Rectangle trackRectangle = new Rectangle(trackX, trackY, trackWidth, trackHeight);
                graphics.setColor(AttributeManager.getInstance().getColor(key, attributeValue));
                graphics.fill(trackRectangle);
                addMousableRegion(new MouseableRegion(trackRectangle, key, attributeValue));
            }
            trackX += trackWidth + AttributeHeaderPanel.COLUMN_BORDER_WIDTH;
        }
        trackY += trackHeight;

        return trackY;
    }*/

    private void init() {

        setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
        if (!PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_ATTRIBUTE_VIEWS_KEY)) {
            setSize(0, getHeight());
        }
        setBackground(new java.awt.Color(255, 255, 255));
        setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
        setPreferredSize(new java.awt.Dimension(0, 0));
        setVerifyInputWhenFocusTarget(false);

        MouseInputAdapter mouseAdapter = new AttributePanelMouseAdapter();
        addMouseMotionListener(mouseAdapter);
        addMouseListener(mouseAdapter);
    }

    /**
     * Method description
     *
     * @param x
     * @param y
     * @return
     */
    public String getPopupMenuTitle(int x, int y) {
        Collection<Track> selectedTracks = getSelectedTracks();
        int selectedTrackCount = selectedTracks.size();
        if (selectedTrackCount == 1) {
            return selectedTracks.iterator().next().getName();
        } else {
            String keyValue = "";
            for (MouseableRegion region : this.getMouseRegions()) {
                if (region.containsPoint(x, y)) {
                    keyValue = region.getText();
                }
            }
            return keyValue + " (" + selectedTrackCount + " tracks)";
        }
    }

    /**
     * Method description
     *
     * @param x
     * @param y
     * @return
     */
    public String getMouseDoc(int x, int y) {

        List<MouseableRegion> mouseRegions = getMouseRegions();

        for (MouseableRegion mr : mouseRegions) {
            if (mr.containsPoint(x, y)) {
                return mr.getText();
            }
        }
        return "";
    }

    /**
     * Method description
     *
     * @return
     */
    public int getAttributeColumnWidth() {
        return AttributeHeaderPanel.ATTRIBUTE_COLUMN_WIDTH;
    }

    // Packable interface


    private int calculatePackWidth() {

        if (!PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_ATTRIBUTE_VIEWS_KEY)) {
            return 0;
        }

        HashSet<String> attributeKeys = new HashSet(AttributeManager.getInstance().getAttributeNames());
        final Set<String> hiddenAttributes = IGV.getInstance().getSession().getHiddenAttributes();
        if (hiddenAttributes != null) attributeKeys.removeAll(hiddenAttributes);

        int attributeCount = attributeKeys.size();
        int packWidth = (attributeCount) * (AttributeHeaderPanel.ATTRIBUTE_COLUMN_WIDTH +
                AttributeHeaderPanel.COLUMN_BORDER_WIDTH) + AttributeHeaderPanel.COLUMN_BORDER_WIDTH;
        return packWidth;
    }

    /**
     * Method description
     *
     * @param x
     * @param y
     * @param width
     * @param height
     */
    @Override
    public void setBounds(int x, int y, int width, int height) {
        super.setBounds(x, y, calculatePackWidth(), height);
    }

    /**
     * Method description
     */
    public void packComponent() {
        int newWidth = calculatePackWidth();

        Dimension dimension = getSize();
        dimension = new Dimension(newWidth, dimension.height);
        setMinimumSize(dimension);
        setMaximumSize(dimension);
        setSize(dimension);
        setPreferredSize(dimension);

    }

    class AttributePanelMouseAdapter extends MouseInputAdapter {

        int lastMousePressX = 0;

        public void mousePressed(final MouseEvent e) {

            if (log.isDebugEnabled()) {
                log.debug("Enter mousePressed");
            }
            clearTrackSelections();
            selectTracks(e);
            if (e.isPopupTrigger()) {
                TrackClickEvent te = new TrackClickEvent(e, null);
                openPopupMenu(te);

            }
            IGV.getInstance().repaintNamePanels();
        }

        public void mouseReleased(MouseEvent e) {
            // Show Popup Menu.  The track selection is cleared afterwards.
            if (e.isPopupTrigger()) {
                TrackClickEvent te = new TrackClickEvent(e, null);
                openPopupMenu(te);
            }

        }

        public void mouseMoved(MouseEvent e) {
            setToolTipText(getMouseDoc(e.getX(), e.getY()));
        }
    }

}
