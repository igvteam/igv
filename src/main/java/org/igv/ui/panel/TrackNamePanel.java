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
import org.igv.track.Track;
import org.igv.track.TrackClickEvent;
import org.igv.track.TrackGroup;
import org.igv.ui.util.IGVMouseInputAdapter;
import org.jdesktop.layout.GroupLayout;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * @author jrobinso
 */
public class TrackNamePanel extends TrackPanelComponent implements Paintable {

    private static Logger log = LogManager.getLogger(TrackNamePanel.class);


    List<GroupExtent> groupExtents = new ArrayList();

    public TrackNamePanel(TrackPanel trackPanel) {
        super(trackPanel);
        init();
    }

    @Override
    public void paintComponent(Graphics g) {

        super.paintComponent(g);

        Rectangle trackRectangle = new Rectangle(getBounds());
        Rectangle visibleRect = getVisibleRect();
        trackRectangle.x = 0; // getBounds returns a rectangle with x/y relative to the parent, we want relative to this component
        trackRectangle.y = 0;

        if (trackRectangle != null && trackRectangle.height > 10) {

            Graphics2D fontGraphics = null;

            try {
                if (darkMode) {
                    setBackground(UIManager.getColor("Panel.background"));
                }

                fontGraphics = (Graphics2D) g.create();

                final Color backgroundColor = darkMode ?
                        UIManager.getColor("Panel.background") :
                        PreferencesManager.getPreferences().getAsColor(Constants.BACKGROUND_COLOR);
                fontGraphics.setBackground(backgroundColor);
                fontGraphics.setColor(backgroundColor);
                fontGraphics.fillRect(visibleRect.x, visibleRect.y, visibleRect.width, visibleRect.height);

                fontGraphics.setColor(darkMode ? Color.white : Color.BLACK);

                if (PreferencesManager.getPreferences().getAntiAliasing()) {
                    ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
                }

                removeMousableRegions();

                paintImpl(fontGraphics, trackRectangle, getVisibleRect(),false);
            } finally {
                fontGraphics.dispose();
            }
        }
    }


    public void paintOffscreen(Graphics2D g, Rectangle rect, boolean batch) {

        Graphics borderGraphics = g.create();
        borderGraphics.setColor(Color.lightGray);

        g.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

        paintImpl(g, rect, rect,true);

        borderGraphics.drawRect(rect.x, rect.y, rect.width - 1, rect.height - 1);
        borderGraphics.dispose();
    }

    @Override
    public int getSnapshotHeight(boolean batch) {
        return getHeight();
    }


    private void paintImpl(Graphics2D g, Rectangle trackRectangle, Rectangle visibleRect, boolean snapshot) {

            Track track = getTrack();
            if (track.isVisible()) {
                track.renderName(g, trackRectangle, visibleRect);
            }
    }

    private void init() {

        //    setBorder(javax.swing.BorderFactory.createLineBorder(Color.black));
        //    setBackground(new java.awt.Color(255, 255, 255));
        GroupLayout dataTrackNamePanelLayout = new GroupLayout(this);
        setLayout(dataTrackNamePanelLayout);
        dataTrackNamePanelLayout.setHorizontalGroup(
                dataTrackNamePanelLayout.createParallelGroup(GroupLayout.LEADING).add(0, 148, Short.MAX_VALUE));
        dataTrackNamePanelLayout.setVerticalGroup(
                dataTrackNamePanelLayout.createParallelGroup(GroupLayout.LEADING).add(0, 528, Short.MAX_VALUE));

        NamePanelMouseAdapter mouseAdapter = new NamePanelMouseAdapter();
        addMouseListener(mouseAdapter);
        addMouseMotionListener(mouseAdapter);

    }


    public String getTooltipTextForLocation(int x, int y) {

        List<MouseableRegion> mouseableRegions = TrackNamePanel.this.getMouseRegions();

        String text = null;
        for (MouseableRegion mouseableRegion : mouseableRegions) {
            if (mouseableRegion.containsPoint(x, y)) {
                Collection<Track> tracks = mouseableRegion.getTracks();
                if (tracks != null && tracks.size() == 1) {
                    Track track = tracks.iterator().next();
                    text = track.getTooltipText(y);
                } else {
                    text = mouseableRegion.getText();
                }
                break;
            }
        }
        return text;
    }



    private TrackGroup getGroup(int y) {
        for (GroupExtent ge : groupExtents) {
            if (ge.contains(y)) {
                return ge.group;
            }
        }
        return null;
    }

    /**
     * Mouse adapter for the track name panel.  Supports multiple selection,
     * popup menu, and drag & drop within or between name panels.
     */
    class NamePanelMouseAdapter extends IGVMouseInputAdapter {

        @Override
        /**
         * Mouse down.  Track selection logic goes here.
         */
        public void mousePressed(MouseEvent e) {
            super.mousePressed(e);

            requestFocus();
            grabFocus();

            if (e.isPopupTrigger()) {
                TrackClickEvent te = new TrackClickEvent(e, null);
                openPopupMenu(te);
            }
        }

        @Override
        public void mouseReleased(MouseEvent e) {

            super.mouseReleased(e);

            if (e.isPopupTrigger()) {
                TrackClickEvent te = new TrackClickEvent(e, null);
                openPopupMenu(te);
            }
        }

        @Override
        public void mouseMoved(MouseEvent e) {
            int x = e.getX();
            int y = e.getY();
            setToolTipText(getTooltipTextForLocation(x, y));
        }

        /**
         * Mouse was clicked.  Delegate action to the track(s) clicked on. .
         *
         * @param e
         */
        @Override
        public void igvMouseClicked(final MouseEvent e) {
            for (MouseableRegion mouseRegion : mouseRegions) {
                if (mouseRegion.containsPoint(e.getX(), e.getY())) {
                    for (Track t : mouseRegion.getTracks()) {
                        t.handleNameClick(e);
                        return;
                    }
                }
            }
        }

    }

    class GroupExtent {
        TrackGroup group;
        int minY;
        int maxY;

        GroupExtent(TrackGroup group, int minY, int maxY) {
            this.group = group;
            this.maxY = maxY;
            this.minY = minY;
        }

        boolean contains(int y) {
            return y > minY && y <= maxY;
        }

        boolean isAfter(int y) {
            return minY > y;
        }
    }

}
