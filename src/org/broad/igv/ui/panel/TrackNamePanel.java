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


import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackClickEvent;
import org.broad.igv.track.TrackGroup;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.dnd.AbstractGhostDropManager;
import org.broad.igv.ui.dnd.GhostDropEvent;
import org.broad.igv.ui.dnd.GhostDropListener;
import org.broad.igv.ui.dnd.GhostGlassPane;
import org.broad.igv.ui.util.UIUtilities;
import org.jdesktop.layout.GroupLayout;

import javax.swing.*;
import javax.swing.event.MouseInputAdapter;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class TrackNamePanel extends TrackPanelComponent implements Paintable {

    private static Logger log = Logger.getLogger(TrackNamePanel.class);


    List<GroupExtent> groupExtents = new ArrayList();
    BufferedImage dndImage = null;
    TrackGroup selectedGroup = null;
    boolean showGroupNames = true;
    boolean showSampleNamesWhenGrouped = false;


    public TrackNamePanel(TrackPanel trackPanel) {
        super(trackPanel);
        init();
    }

    Collection<TrackGroup> getGroups() {
        return getTrackPanel().getGroups();
    }

    private boolean isGrouped() {
        return getGroups().size() > 1;
    }


    @Override
    public void paintComponent(Graphics g) {

        super.paintComponent(g);
        removeMousableRegions();
        Rectangle visibleRect = getVisibleRect();
        paintImpl(g, visibleRect);
    }


    public void paintOffscreen(Graphics2D g, Rectangle rect) {
        g.setColor(Color.white);
        g.fill(rect);
        paintImpl(g, rect);

        Color c = g.getColor();
        g.setColor(Color.darkGray);
        g.drawRect(rect.x, rect.y, rect.width, rect.height);
        g.setColor(c);            //super.paintBorder(g);
        //super.paintBorder(g);
    }


    private void paintImpl(Graphics g, Rectangle visibleRect) {
        // Get available tracks
        Collection<TrackGroup> groups = getGroups();
        boolean isGrouped = groups.size() > 1;


        if (!groups.isEmpty()) {
            final Graphics2D graphics2D = (Graphics2D) g.create();
            graphics2D.setColor(Color.BLACK);
            if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.ENABLE_ANTIALISING)) {
                graphics2D.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
            }

            final Graphics2D greyGraphics = (Graphics2D) g.create();
            greyGraphics.setColor(UIConstants.LIGHT_GREY);

            int regionY = 0;

            groupExtents.clear();

            //Rectangle clipRect = g.getClipBounds();

            for (Iterator<TrackGroup> groupIter = groups.iterator(); groupIter.hasNext(); ) {
                TrackGroup group = groupIter.next();

                if (regionY > visibleRect.getMaxY()) {
                    break;
                }

                if (group.isVisible()) {
                    if (isGrouped) {
                        if (regionY + UIConstants.groupGap >= visibleRect.y && regionY < visibleRect.getMaxY()) {
                            greyGraphics.fillRect(0, regionY + 1, getWidth(), UIConstants.groupGap - 1);
                        }
                        regionY += UIConstants.groupGap;
                    }

                    if (group.isDrawBorder() && regionY + UIConstants.groupGap >= visibleRect.y &&
                            regionY < visibleRect.getMaxY()) {
                        g.drawLine(0, regionY - 1, getWidth(), regionY - 1);
                    }

                    int h = group.getHeight();
                    Rectangle groupRect = new Rectangle(visibleRect.x, regionY, visibleRect.width, h);
                    Rectangle displayableRect = getDisplayableRect(groupRect, visibleRect);
                    regionY = printTrackNames(group, displayableRect, visibleRect, graphics2D, 0, regionY);

                    if (isGrouped) {
                        groupExtents.add(new GroupExtent(group, groupRect.y, groupRect.y + groupRect.height));
                        if (showGroupNames) {
                            //Rectangle displayableRect = getDisplayableRect(groupRect, visibleRect);
                            group.renderName(graphics2D, displayableRect, group == selectedGroup);
                        }
                    }

                    if (group.isDrawBorder()) {
                        g.drawLine(0, regionY, getWidth(), regionY);
                    }
                }

            }
        }
    }

    private Rectangle getDisplayableRect(Rectangle trackRectangle, Rectangle visibleRect) {
        Rectangle rect = null;
        if (visibleRect != null) {
            Rectangle intersectedRect = trackRectangle.intersection(visibleRect);
            if (intersectedRect.height > 15) {
                rect = intersectedRect;
            } else {
                rect = new Rectangle(trackRectangle);
            }
        }
        return rect;

    }

    private int printTrackNames(TrackGroup group, Rectangle visibleRect, Rectangle clipRect,
                                Graphics2D graphics2D, int regionX, int regionY) {


        List<Track> tmp = new ArrayList(group.getTracks());
        final Color backgroundColor = PreferenceManager.getInstance().getAsColor(PreferenceManager.BACKGROUND_COLOR);
        graphics2D.setBackground(backgroundColor);
        graphics2D.clearRect(visibleRect.x, visibleRect.y, visibleRect.width, visibleRect.height);

        for (Track track : tmp) {
            if (track == null) continue;
            track.setY(regionY);
            int trackHeight = track.getHeight();
            if (track.isVisible()) {

                if (regionY + trackHeight >= clipRect.y && regionY < clipRect.getMaxY()) {
                    int width = getWidth();
                    int height = track.getHeight();

                    Rectangle region = new Rectangle(regionX, regionY, width, height);
                    addMousableRegion(new MouseableRegion(region, track));

                    if (!isGrouped() || showSampleNamesWhenGrouped) {
                        Rectangle rect = new Rectangle(regionX, regionY, width, height);
                        //Graphics2D g2D = graphics; //(Graphics2D) graphics.create();
                        if (track.isSelected()) {
                            graphics2D.setBackground(Color.LIGHT_GRAY);
                            graphics2D.clearRect(rect.x, rect.y, rect.width, rect.height);
                        } else {
                            graphics2D.setBackground(backgroundColor);
                        }
                        track.renderName(graphics2D, rect, visibleRect);
                    }

                }
                regionY += trackHeight;
            }
        }
        return regionY;
    }


    private void init() {

        setBorder(javax.swing.BorderFactory.createLineBorder(Color.black));
        setBackground(new java.awt.Color(255, 255, 255));
        GroupLayout dataTrackNamePanelLayout = new org.jdesktop.layout.GroupLayout(this);
        setLayout(dataTrackNamePanelLayout);
        dataTrackNamePanelLayout.setHorizontalGroup(
                dataTrackNamePanelLayout.createParallelGroup(GroupLayout.LEADING).add(0, 148, Short.MAX_VALUE));
        dataTrackNamePanelLayout.setVerticalGroup(
                dataTrackNamePanelLayout.createParallelGroup(GroupLayout.LEADING).add(0, 528, Short.MAX_VALUE));

        NamePanelMouseAdapter mouseAdapter = new NamePanelMouseAdapter();
        addMouseListener(mouseAdapter);
        addMouseMotionListener(mouseAdapter);

        DropListener dndListener = new DropListener(this);
        addGhostDropListener(dndListener);
    }


    @Override
    protected void openPopupMenu(TrackClickEvent te) {

        ArrayList<Component> extraItems = null;
        if (isGrouped()) {
            extraItems = new ArrayList();

            final JMenuItem item = new JCheckBoxMenuItem("Show group names");
            item.setSelected(showGroupNames);
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    showGroupNames = item.isSelected();
                    repaint();
                }
            });
            extraItems.add(item);

            final JMenuItem item2 = new JCheckBoxMenuItem("Show sample names");
            item2.setSelected(showSampleNamesWhenGrouped);
            item2.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    showSampleNamesWhenGrouped = item2.isSelected();
                    repaint();
                }
            });
            extraItems.add(item2);
        }

        super.openPopupMenu(te, extraItems);


    }


    public String getTooltipTextForLocation(int x, int y) {

        List<MouseableRegion> mouseableRegions = TrackNamePanel.this.getMouseRegions();

        String text = null;
        for (MouseableRegion mouseableRegion : mouseableRegions) {
            if (mouseableRegion.containsPoint(x, y)) {
                Collection<Track> tracks = mouseableRegion.getTracks();
                if (tracks != null && tracks.size() == 1) {
                    Track track = tracks.iterator().next();
                    text = track.getNameValueString(y);
                } else {
                    text = mouseableRegion.getText();
                }
                break;
            }
        }
        return text;
    }


    private synchronized void createDnDImage() {
        dndImage = new BufferedImage(getWidth(), 2, BufferedImage.TYPE_INT_ARGB);
        Graphics g = dndImage.getGraphics();
        g.setColor(Color.blue);
        g.drawLine(1, 0, getWidth() - 2, 0);
        g.drawLine(1, 1, getWidth() - 2, 1);

    }

    /**
     * Shift-click,  used to select a range of tracks.
     *
     * @param e
     */
    protected void shiftSelectTracks(MouseEvent e) {
        for (MouseableRegion mouseRegion : mouseRegions) {
            if (mouseRegion.containsPoint(e.getX(), e.getY())) {
                Collection<Track> clickedTracks = mouseRegion.getTracks();
                if (clickedTracks != null && clickedTracks.size() > 0) {
                    Track t = clickedTracks.iterator().next();
                    IGV.getInstance().shiftSelectTracks(t);
                }
                return;
            }
        }
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
    class NamePanelMouseAdapter extends MouseInputAdapter {

        boolean isDragging = false;
        List<Track> dragTracks = new ArrayList();
        Point dragStart = null;

        @Override
        /**
         * Mouse down.  Track selection logic goes here.
         */
        public void mousePressed(MouseEvent e) {

            dragStart = e.getPoint();

            requestFocus();
            grabFocus();

            boolean isGrouped = isGrouped();

            if (e.isPopupTrigger()) {
                if (isGrouped) {
                    clearTrackSelections();
                    TrackGroup g = getGroup(e.getY());
                    if (null == g || g == selectedGroup) {
                        selectedGroup = null;
                    } else {
                        selectGroup(g);
                    }
                } else if (!isTrackSelected(e)) {
                    clearTrackSelections();
                    selectTracks(e);
                }
                TrackClickEvent te = new TrackClickEvent(e, null);
                openPopupMenu(te);
            } // meta (mac) or control,  toggle selection]
            else if (e.getButton() == MouseEvent.BUTTON1) {

                if (isGrouped) {
                    clearTrackSelections();
                    TrackGroup g = getGroup(e.getY());
                    if (g == selectedGroup) {
                        selectedGroup = null;
                    } else {
                        selectGroup(getGroup(e.getY()));
                    }
                } else {
                    if (e.isMetaDown() || e.isControlDown()) {
                        toggleTrackSelections(e);
                    } else if (e.isShiftDown()) {
                        shiftSelectTracks(e);
                    } else if (!isTrackSelected(e)) {
                        clearTrackSelections();
                        selectTracks(e);
                    }
                }
            } else {
                if (isGrouped) {

                } else if (!isTrackSelected(e)) {
                    clearTrackSelections();
                    selectTracks(e);
                }
            }


            IGV.getInstance().repaintNamePanels();

        }

        public void mouseReleased(MouseEvent e) {

            if (log.isTraceEnabled()) {
                log.trace("Enter mouseReleased");
            }

            if (isDragging) {


                Component c = e.getComponent();

                IGV.getInstance().endDnD();
                GhostGlassPane glassPane = IGV.getInstance().getDnDGlassPane();

                Point p = (Point) e.getPoint().clone();
                SwingUtilities.convertPointToScreen(p, c);

                Point eventPoint = (Point) p.clone();
                SwingUtilities.convertPointFromScreen(p, glassPane);

                glassPane.setPoint(p);
                glassPane.setVisible(false);
                glassPane.setImage(null);

                fireGhostDropEvent(new GhostDropEvent(dragStart, eventPoint, dragTracks));

                if (selectedGroup != null) {
                    int idx = getGroupGapNumber(e.getY());
                    TrackPanel dataTrackView = (TrackPanel) getParent();
                    dataTrackView.moveGroup(selectedGroup, idx);
                    dataTrackView.repaint();
                }
                selectedGroup = null;


            }

            if (e.isPopupTrigger()) {
                TrackClickEvent te = new TrackClickEvent(e, null);
                openPopupMenu(te);
            } else {
                if (!isDragging && !e.isMetaDown() && !e.isControlDown() &&
                        !e.isShiftDown()) {
                    clearTrackSelections();
                    selectTracks(e);
                    IGV.getInstance().repaintNamePanels();
                }
            }

            isDragging = false;
            dragTracks.clear();
            dndImage = null;


        }


        public void mouseDragged(MouseEvent e) {

            Component c = e.getComponent();
            if (e.isPopupTrigger()) {
                return;
            }
            if (!isDragging) {

                if (dragStart == null) {
                    dragStart = e.getPoint();
                    return;
                } else if (e.getPoint().distance(dragStart) < 5) {
                    return;
                }

                dragStart.x = getWidth() / 2;
                IGV.getInstance().startDnD();

                if (dndImage == null) {
                    createDnDImage();
                }
                IGV.getInstance().getDnDGlassPane().setImage(dndImage);
                isDragging = true;
                dragTracks.clear();
                dragTracks.addAll(IGV.getInstance().getSelectedTracks());


                if (getGroups().size() > 0) {
                    selectedGroup = getGroup(e.getY());
                } else {
                    selectedGroup = null;
                }

                // Code below paints target component on the dndImage.  It needs modified to paint some representation
                // of the selectect tracks, probably the track names printed as a list.
            }
            if (isDragging) {

                final GhostGlassPane glassPane = IGV.getInstance().getDnDGlassPane();

                Point p = (Point) e.getPoint().clone();
                p.x = getWidth() / 2;
                SwingUtilities.convertPointToScreen(p, c);
                SwingUtilities.convertPointFromScreen(p, glassPane);

                glassPane.setPoint(p);

                UIUtilities.invokeOnEventThread(new Runnable() {

                    public void run() {
                        Rectangle bounds = new Rectangle(getBounds());
                        bounds.height = 10000;
                        glassPane.paintImmediately(bounds);
                    }
                });
            }
        }

        @Override
        public void mouseMoved(MouseEvent e) {
            int x = e.getX();
            int y = e.getY();
            setToolTipText(getTooltipTextForLocation(x, y));
        }


        /**
         * Mouse was clicked.  Delegate single-click action to the track(s) clicked on.   We won't know if this
         * is a double click or not until the double-click interval has passed, so defer the action with a
         * TimerTask.  If a second click arrives it will be canceled.
         *
         * @param e
         */
        @Override
        public void mouseClicked(final MouseEvent e) {

            // If this is the second click of a double click, cancel the scheduled single click task.
            if (e.getClickCount() > 1) {
                clickScheduler.cancelClickTask();
                return;
            }

            TimerTask clickTask = new TimerTask() {

                @Override
                public void run() {
                    for (MouseableRegion mouseRegion : mouseRegions) {
                        if (mouseRegion.containsPoint(e.getX(), e.getY())) {
                            for (Track t : mouseRegion.getTracks()) {
                                t.handleNameClick(e);
                            }
                            return;
                        }
                    }//To change body of implemented methods use File | Settings | File Templates.
                }
            };
            //clickScheduler.scheduleClickTask(clickTask);
            clickTask.run();
        }

        protected void fireGhostDropEvent(GhostDropEvent evt) {
            Iterator it = TrackNamePanel.dropListeners.iterator();
            while (it.hasNext()) {
                ((GhostDropListener) it.next()).ghostDropped(evt);
            }
        }
    }


    class DropListener extends AbstractGhostDropManager {

        TrackNamePanel panel;

        public DropListener(TrackNamePanel target) {
            super(target);
            this.panel = target;

        }

        public void ghostDropped(GhostDropEvent e) {
            Point startPoint = e.getStartLocation();
            Point dropPoint = getTranslatedPoint(e.getDropLocation());


            Rectangle bounds = component.getVisibleRect();
            boolean isInTarget = dropPoint.y > bounds.y && dropPoint.y < bounds.getMaxY();

            if (isInTarget) {
                tracksDropped(startPoint, dropPoint, e.getTracks());
                e.removeTracksFromSource();
                e.setTracksDropped(true);
            } else {
                TrackPanel view = ((TrackPanel) getParent());
                if (e.isTracksDropped()) {
                    view.removeTracks(e.getTracks());
                } else {
                    // Defer removal until we are sure the tracks are dropped in another panel
                    e.addSourcePanel(view);
                }
            }
        }

        void tracksDropped(Point startPoint, Point dropPoint, List<Track> tracks) {

            // This cast is horrid but we can't fix everything at once.
            TrackPanel panel = ((TrackPanel) getParent());
            List<MouseableRegion> regions = getMouseRegions();

            if (regions.isEmpty()) {
                // empty panel,  just add the tracks
                panel.addTracks(tracks);
            } else {
                // Find the regions containing the startPoint and point
                boolean before = true;
                MouseableRegion dropRegion = null;
                MouseableRegion startRegion = null;
                for (MouseableRegion region : regions) {
                    if (region.containsPoint(dropPoint.x, dropPoint.y)) {
                        dropRegion = region;
                        Rectangle bnds = dropRegion.getBounds();
                        int dy1 = (dropPoint.y - bnds.y);
                        int dy2 = bnds.height - dy1;
                        before = dy1 < dy2;
                    }
                    if (region.containsPoint(startPoint.x, startPoint.y)) {
                        startRegion = region;
                    }
                    if (dropRegion != null && startRegion != null) {
                        break;
                    }
                }


                Track dropTrack = null;
                if (dropRegion != null) {
                    Iterator<Track> tmp = dropRegion.getTracks().iterator();
                    if (tmp.hasNext()) {
                        dropTrack = tmp.next();
                    }
                }
                panel.moveSelectedTracksTo(tracks, dropTrack, before);
            }


        }
    }

    private void selectGroup(TrackGroup group) {
        selectedGroup = group;
        if (selectedGroup != null) {
            for (Track t : selectedGroup.getTracks()) {
                t.setSelected(true);
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

    int getGroupGapNumber(int y) {
        for (int i = 0; i < groupExtents.size(); i++) {
            if (groupExtents.get(i).isAfter(y)) {
                return i;
            }
        }
        return groupExtents.size();
    }


    // Track D&D support follows
    // TODO -- this use of a static is really bad,  bugs and memory leaks waiting to happen.  Redesign this.
    static List<DropListener> dropListeners = new ArrayList();


    private static void addGhostDropListener(DropListener listener) {
        if (listener != null) {
            dropListeners.add(listener);
        }
    }

    public static void removeDropListenerFor(TrackNamePanel panel) {
        List<DropListener> removeThese = new ArrayList();
        for (DropListener dl : dropListeners) {
            if (dl.panel == panel) {
                removeThese.add(dl);
            }
        }
        dropListeners.removeAll(removeThese);
    }


}
