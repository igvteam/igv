/*
 * TrackPanel.java
 *
 * Created on Sep 5, 2007, 4:09:39 PM
 *
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.igv.ui.panel;


import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.track.AttributeManager;
import org.igv.track.Track;
import org.igv.track.TrackClickEvent;
import org.igv.ui.IGV;

import javax.swing.event.MouseInputAdapter;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.List;
import java.util.Set;

/**
 * @author jrobinso
 */
public class AttributePanel extends TrackPanelComponent implements Paintable {

    private static Logger log = LogManager.getLogger(AttributePanel.class);

    /**
     * Constructs ...
     */
    public AttributePanel(TrackPanel trackPanel) {
        super(trackPanel);
        init();
    }

    private void init() {

        if (!PreferencesManager.getPreferences().getAsBoolean(Constants.SHOW_ATTRIBUTE_VIEWS_KEY)) {
            setSize(0, getHeight());
        }
        setBorder(null);
        setPreferredSize(new java.awt.Dimension(0, 0));
        setVerifyInputWhenFocusTarget(false);

        MouseInputAdapter mouseAdapter = new AttributePanelMouseAdapter();
        addMouseMotionListener(mouseAdapter);
        addMouseListener(mouseAdapter);
    }

    @Override
    protected void paintComponent(Graphics g) {

        super.paintComponent(g);
        Rectangle trackRectangle = new Rectangle(getBounds());
        trackRectangle.x = 0;
        trackRectangle.y = 0;
        removeMousableRegions();
        paintImpl((Graphics2D) g, trackRectangle);

    }

    @Override
    public void paintOffscreen(Graphics2D g, Rectangle rect, boolean batch) {
        paintImpl(g, rect);
    }

    public void paintImpl(Graphics2D g, Rectangle rect) {

        Track track = getTrack();
        if (track.isVisible()) {
            List<String> names = AttributeManager.getInstance().getAttributeNames();

            if (names == null) {
                return;
            }

            final Set<String> hiddenAttributes = IGV.getInstance().getSession().getHiddenAttributes();
            if (hiddenAttributes != null) names.removeAll(hiddenAttributes);

            if (names.size() > 0 && track.getSampleGroups() != null) {
                track.renderAttributes(g, rect, names, mouseRegions);
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
        Track track = getTrack();
        if (track != null) {
            return track.getName();
        }
        String keyValue = "";
        for (MouseableRegion region : this.getMouseRegions()) {
            if (region.containsPoint(x, y)) {
                keyValue = region.getText();
            }
        }
        return keyValue;
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
            if (e.isPopupTrigger()) {
                TrackClickEvent te = new TrackClickEvent(e, null);
                openPopupMenu(te);

            }
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

        @Override
        public void mouseEntered(MouseEvent e) {
            setToolTipText(getMouseDoc(e.getX(), e.getY()));
        }
    }

}
