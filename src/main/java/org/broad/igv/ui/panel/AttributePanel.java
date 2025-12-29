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

import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
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
import java.util.List;
import java.util.*;

/**
 * @author jrobinso
 */
public class AttributePanel extends TrackPanelComponent implements  Paintable {

    private static Logger log = LogManager.getLogger(AttributePanel.class);


    /**
     * Constructs ...
     */
    public AttributePanel(TrackPanel trackPanel) {
        super(trackPanel);
        init();
    }

    private void init() {

        setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
        if (!PreferencesManager.getPreferences().getAsBoolean(Constants.SHOW_ATTRIBUTE_VIEWS_KEY)) {
            setSize(0, getHeight());
        }
        setBackground(new java.awt.Color(255, 255, 255));
        setBorder(javax.swing.BorderFactory.createLineBorder(Color.lightGray));
        setPreferredSize(new java.awt.Dimension(0, 0));
        setVerifyInputWhenFocusTarget(false);

        MouseInputAdapter mouseAdapter = new AttributePanelMouseAdapter();
        addMouseMotionListener(mouseAdapter);
        addMouseListener(mouseAdapter);
    }

    @Override
    protected void paintComponent(Graphics g) {

        super.paintComponent(g);
        Rectangle visibleRect = getVisibleRect();
        removeMousableRegions();
        paintImpl((Graphics2D) g, visibleRect, false);

    }

    @Override
    public void paintOffscreen(Graphics2D g, Rectangle rect, boolean batch) {
        paintImpl(g, rect, batch);
        Graphics2D borderGraphics = (Graphics2D) g.create();
        borderGraphics.setColor(Color.lightGray);
        rect.height -=1;
        borderGraphics.draw(rect);
        borderGraphics.dispose();
    }

    public void paintImpl(Graphics2D g, Rectangle rect, boolean batch) {

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

                        for (Track track : group.getVisibleTracks()) {
                            if (track == null) continue;
                            int trackHeight = track.getHeight();
                            if (regionY > bottom) {
                                break;
                            }

                            if (track.isVisible()) {
                                int border = trackHeight < 5 ? 0 : 1;
                                if (regionY + trackHeight >= rect.y) {
                                    Rectangle trackRectangle = new Rectangle(left, regionY + border, getWidth(), trackHeight - border);
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

                // Border between columns
                final Graphics2D columnBorderGraphics = (Graphics2D) g.create();
                columnBorderGraphics.setColor(Color.lightGray);
                final int colWidth = AttributeHeaderPanel.ATTRIBUTE_COLUMN_WIDTH + AttributeHeaderPanel.COLUMN_BORDER_WIDTH;
                for (int x = 1; x < rect.x + rect.width; x += colWidth) {
                    columnBorderGraphics.fillRect(
                            x + AttributeHeaderPanel.ATTRIBUTE_COLUMN_WIDTH, rect.y, AttributeHeaderPanel.COLUMN_BORDER_WIDTH, rect.height);
                }

                graphics2D.dispose();
                greyGraphics.dispose();
                columnBorderGraphics.dispose();

            }

        }

    }

    @Override
    public int getSnapshotHeight(boolean batch) {
        return getHeight();
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
