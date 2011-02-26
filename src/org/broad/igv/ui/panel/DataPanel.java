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

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackClickEvent;
import org.broad.igv.track.TrackGroup;
import org.broad.igv.ui.*;
import org.broad.igv.ui.util.DataPanelTool;

import javax.swing.*;
import javax.swing.event.MouseInputAdapter;
import java.awt.*;
import java.awt.event.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.List;

/**
 * The main panel for displaying tracks and data.  A DataPanel is always associated with a ReferenceFrame.  Normally
 * there is a single reference frame (and thus panel), but when "gene list" or other split screen views are
 * invoked there can be multiple panels.
 *
 * @author jrobinso
 */
public class DataPanel extends JComponent implements Paintable {

    private static Logger log = Logger.getLogger(DataPanel.class);

    // TODO move this to some central place
    final static private boolean IS_MAC = System.getProperty("os.name").toLowerCase().startsWith("mac");

    private DataPanelTool defaultTool;
    private DataPanelTool currentTool;
    private Point tooltipTextPosition;
    private ReferenceFrame frame;
    DataPanelContainer parent;
    private DataPanelPainter painter;


    public DataPanel(ReferenceFrame frame, DataPanelContainer parent) {
        init();
        this.defaultTool = new PanAndZoomTool(this);
        this.currentTool = defaultTool;
        this.frame = frame;
        this.parent = parent;
        this.setBackground(Color.white);
        setFocusable(true);
        setAutoscrolls(true);
        setToolTipText("");
        painter = new DataPanelPainter();
        //add a listener that kills the tooltip when we navigate away from tne window
        addFocusListener(new DataPanelFocusListener(this));
    }

    /**
     * This class exists entirely to tell the ToolTip to go away when the DataPanel loses focus.
     * Without this, the tooltip stays up and annoyingly visible even if another window is on top of IGV
     * <p/>
     * It's also necessary to let the tooltip manager know when we get focus back.  Otherwise, the
     * tooltip won't display until you move the mouse out of the panel and back.
     */
    private class DataPanelFocusListener implements FocusListener {
        protected Component cmp;

        public DataPanelFocusListener(Component cmp) {
            this.cmp = cmp;
        }

        public void focusGained(FocusEvent focusEvent) {
            //This is a bit of a hack -- if the mouse is in the panel, generate a mouseEntered event to tell
            // the tooltipmanager that the mouse is back, since we told it earlier that it went away when it hadn't
            if (getMousePosition() != null && cmp.contains(getMousePosition()))
                ToolTipManager.sharedInstance().mouseEntered(new MouseEvent(cmp, 0, 0, 0,
                        (int) getMousePosition().getX(), (int) getMousePosition().getY(), 0, false));

        }

        public void focusLost(FocusEvent focusEvent) {
            //This is a bit hacky, but it works and shouldn't cause any harm
            ToolTipManager.sharedInstance().mouseExited(new MouseEvent(cmp, 0, 0, 0, 0, 0, 0, true));
        }
    }


    /**
     * @return
     */
    public JScrollBar getVerticalScrollbar() {
        Component sp = getParent();
        while (sp != null && !(sp instanceof JScrollPane)) {
            sp = sp.getParent();
        }
        return sp == null ? null : ((JScrollPane) sp).getVerticalScrollBar();
    }


    public void setCurrentTool(final AbstractDataPanelTool tool) {
        this.currentTool = (tool == null) ? defaultTool : tool;
        if (currentTool != null) {
            setCursor(currentTool.getCursor());
        }
    }


    @Override
    public void paintComponent(final Graphics g) {


        super.paintComponent(g);

        RenderContext context = null;
        try {

            Rectangle clipBounds = g.getClipBounds();
            final Rectangle damageRect = clipBounds == null ? getVisibleRect() : clipBounds.intersection(getVisibleRect());
            Graphics2D graphics2D = (Graphics2D) g; //(Graphics2D) g.create();
            String genomeId = GenomeManager.getInstance().getGenomeId();

            context = new RenderContext(genomeId, this, graphics2D, frame, this.getVisibleRect());

            if (IS_MAC) {
                this.applyMacPerformanceHints((Graphics2D) g);
            }

            graphics2D.setBackground(getBackground());
            graphics2D.clearRect(damageRect.x, damageRect.y, damageRect.width, damageRect.height);

            final Collection<TrackGroup> groups = new ArrayList(parent.getTrackGroups());

            // If there are no tracks to paint, just exit
            boolean hasTracks = false;
            for (TrackGroup group : groups) {
                if (group.getTracks().size() > 0) {
                    hasTracks = true;
                    break;
                }
            }
            if (!hasTracks) {
                removeMousableRegions();
                return;
            }


            int trackWidth = getWidth();
            int trackHeight = getHeight();
            painter.paint(groups, context, trackWidth, trackHeight, getBackground(), damageRect);

            // Compute mouse senstive regions
            computeMousableRegions(groups);


            // If there is a partial ROI in progress draw it first
            if (currentTool instanceof RegionOfInterestTool) {
                int startLoc = ((RegionOfInterestTool) currentTool).getRoiStart();
                if (startLoc > 0) {
                    int start = frame.getScreenPosition(startLoc);
                    g.setColor(Color.BLACK);
                    graphics2D.drawLine(start, 0, start, getHeight());
                }
            }

            if (IGVMainFrame.getInstance().isShowRegionsOfInterestBarsOn()) {
                drawAllRegions(g);
            }

        } finally {

            if (context != null) {
                context.dispose();
            }
        }
    }

    /**
     * Paint method designed to paint to an offscreen image
     *
     * @param g
     * @param rect
     */

    public void paintOffscreen(final Graphics2D g, Rectangle rect) {

        RenderContext context = null;
        try {

            String genomeId = GenomeManager.getInstance().getGenomeId();

            context = new RenderContext(genomeId, null, g, frame, rect);
            final Collection<TrackGroup> groups = new ArrayList(parent.getTrackGroups());
            int width = rect.width;
            int height = rect.height;
            painter.paint(groups, context, width, height, getBackground(), rect);

            super.paintBorder(g);

        } finally {

            if (context != null) {
                context.dispose();
            }
        }
    }


    /**
     * Draw vertical lines demarcating regions of interest.
     */
    public void drawAllRegions(final Graphics g) {

        // TODO -- get rid of this ugly reference to IGVMainFrame
        Collection<RegionOfInterest> regions =
                IGVMainFrame.getInstance().getSession().getRegionsOfInterest(frame.getChrName());

        if ((regions == null) || regions.isEmpty()) {
            return;
        }

        boolean drawBars = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_REGION_BARS);
        Graphics2D graphics2D = (Graphics2D) g.create();
        try {

            for (RegionOfInterest regionOfInterest : regions) {
                if (drawBars || regionOfInterest == RegionOfInterestPanel.getSelectedRegion()) {
                    drawRegion(graphics2D, regionOfInterest);
                }
            }
        } finally {
            if (graphics2D != null) {
                graphics2D.dispose();
            }
        }
    }

    private boolean drawRegion(Graphics2D graphics2D, RegionOfInterest regionOfInterest) {
        Integer regionStart = regionOfInterest.getStart();
        if (regionStart == null) {
            return true;
        }

        Integer regionEnd = regionOfInterest.getEnd();
        if (regionEnd == null) {
            regionEnd = regionStart;
        }
        ReferenceFrame referenceFrame = frame;
        int start = referenceFrame.getScreenPosition(regionStart);
        int end = referenceFrame.getScreenPosition(regionEnd);

        // Set foreground color of boundaries
        int height = getHeight();
        graphics2D.setColor(regionOfInterest.getForegroundColor());
        graphics2D.drawLine(start, 0, start, height);
        graphics2D.drawLine(end, 0, end, height);
        return false;
    }

    protected String generateTileKey(final String chr, int t,
                                     final int zoomLevel) {

        // Fetch image for this chromosome, zoomlevel, and tile.  If found
        // draw immediately
        final String key = chr + "_z_" + zoomLevel + "_t_" + t;
        return key;
    }

    /**
     * TODO -- it is assumed the tracks are drawn in the same order for the
     * image.  This is brittle, the rects could be computed when the image
     * is drawn.  Really the y position and height are all that is needed.
     *
     * @param groups
     */
    private void computeMousableRegions(Collection<TrackGroup> groups) {

        removeMousableRegions();

        int trackX = 0;
        int trackY = 0;

        for (Iterator<TrackGroup> groupIter = groups.iterator(); groupIter.hasNext();) {
            TrackGroup group = groupIter.next();

            if (groups.size() > 1) {
                trackY += UIConstants.groupGap;
            }
            for (Track track : group.getTracks()) {
                if (track == null) continue;
                if (track.isVisible()) {

                    int trackWidth = getWidth();
                    int trackHeight = track.getHeight();

                    Rectangle actualAreaOnDataPanel = new Rectangle(trackX,
                            trackY, (trackWidth),
                            (trackHeight));

                    addMousableRegion(new MouseableRegion(actualAreaOnDataPanel, track));

                    trackY += track.getHeight();
                }
            }
        }

    }


    /**
     * Do not remove - Used for debugging only
     *
     * @param trackName
     */
    public void debugDump(String trackName) {

        // Get the view that holds the track name, attribute and data panels
        TrackPanel trackView = (TrackPanel) getParent();

        if (trackView == null) {
            return;
        }


        if (trackView.hasTracks()) {
            String name = parent.getTrackSetID().toString();
            System.out.println(
                    "\n\n" + name + " Track COUNT:" + trackView.getTracks().size());
            System.out.println(
                    "\t\t\t\t" + name + " scrollpane height     = " + trackView.getScrollPane().getHeight());
            System.out.println(
                    "\t\t\t\t" + name + " viewport height       = " + trackView.getViewportHeight());
            System.out.println(
                    "\t\t\t\t" + name + " TrackView min height  = " + trackView.getMinimumSize().getHeight());
            System.out.println(
                    "\t\t\t\t" + name + " TrackView pref height = " + trackView.getPreferredSize().getHeight());
            System.out.println(
                    "\t\t\t\t" + name + " TrackView height      = " + trackView.getSize().getHeight());
        }

    }

    /**
     * Return html formatted text for mouse position (pixels).
     * TODO  this will be a lot easier when each track has its own panel.
     */
    static DecimalFormat locationFormatter = new DecimalFormat();

    /**
     * Method description
     *
     * @param x
     * @param y
     * @return
     */
    public Track getTrack(int x, int y) {
        for (MouseableRegion mouseRegion : parent.getTrackRegions()) {
            if (mouseRegion.containsPoint(x, y)) {
                return mouseRegion.getTracks().iterator().next();
            }
        }
        return null;

    }

    /**
     * Method description
     *
     * @param x
     * @param y
     */
    public void updateTooltipText(int x, int y) {
        double location = frame.getChromosomePosition(x);
        double displayLocation = location + 1;

        Track track = null;
        List<MouseableRegion> regions = parent.getTrackRegions();
        StringBuffer popupTextBuffer = new StringBuffer();
        for (MouseableRegion mouseRegion : regions) {
            if (mouseRegion.containsPoint(x, y)) {
                track = mouseRegion.getTracks().iterator().next();
                if (track != null) {

                    // Calculate y relative to track bounds
                    //int yRelative = y - mouseRegion.getBounds().y;

                    // First see if there is an overlay track.  If there is, give
                    // it first crack
                    List<Track> overlays = IGVMainFrame.getInstance().getTrackManager().getOverlayTracks(track);

                    if (overlays != null) {
                        for (Track overlay : overlays) {
                            if ((overlay != track) && (overlay.getValueStringAt(
                                    frame.getChrName(), displayLocation, y, frame) != null)) {
                                popupTextBuffer.append(getPopUpText(overlay, displayLocation, y));
                                popupTextBuffer.append("<br>");
                            }
                        }
                    }

                    popupTextBuffer.append(getPopUpText(track, displayLocation, y));
                }
                break;
            }
        }

        String puText = popupTextBuffer.toString().trim();
        if (puText.length() > 0) {

            //popupTextBuffer.append("<br>Location: ");
            //popupTextBuffer.append(frame.getChrName() + ":" + locationFormatter.format((int) displayLocation));

            setToolTipText("<html>" + puText);
        }
    }

    private boolean isWaitingForToolTipText = false;


    @Override
    final public String getToolTipText() {

        if (!isWaitingForToolTipText) {
            isWaitingForToolTipText = true;
            if (tooltipTextPosition != null) {
                updateTooltipText(tooltipTextPosition.x, tooltipTextPosition.y);
            }
            isWaitingForToolTipText = false;
        }
        return super.getToolTipText();
    }


    public String getPopupMenuTitle(int x, int y) {

        Collection<Track> tracks = parent.getSelectedTracks();
        String popupTitle = "";

        // Title for the popup
        if (!tracks.isEmpty()) {
            for (Track track : tracks) {
                popupTitle = track.getName();
                break;
            }
        }
        return popupTitle;
    }

    /**
     * @param track
     * @param location in genomic coordinates
     * @param y        pixel position in panel coordinates
     * @return
     */
    String getPopUpText(Track track, double location, int y) {

        StringBuffer buf = new StringBuffer();
        String value = track.getValueStringAt(frame.getChrName(), location, y, frame);
        if (value != null) {
            buf.append(value);
        }
        return buf.toString();
    }


    private void init() {
        setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
        setBackground(new java.awt.Color(255, 255, 255));
        setRequestFocusEnabled(false);

        // Key Events
        KeyAdapter keyAdapter = new KeyAdapter() {

            @Override
            public void keyPressed(KeyEvent e) {
                if (e.getKeyChar() == '+' || e.getKeyCode() == KeyEvent.VK_PLUS) {
                    WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();
                    try {
                        frame.incrementZoom(1);
                    } finally {
                        WaitCursorManager.removeWaitCursor(token);
                    }
                } else if (e.getKeyChar() == '-' || e.getKeyCode() == KeyEvent.VK_PLUS) {
                    WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();
                    try {
                        frame.incrementZoom(-1);
                    } finally {
                        WaitCursorManager.removeWaitCursor(token);
                    }

                } else if (e.getKeyCode() == KeyEvent.VK_RIGHT) {
                    frame.shiftOriginPixels(5);
                } else if (e.getKeyCode() == KeyEvent.VK_LEFT) {
                    frame.shiftOriginPixels(-5);
                } else if (e.getKeyCode() == KeyEvent.VK_HOME) {
                    WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();
                    try {
                        frame.shiftOriginPixels(-getWidth());
                        frame.recordHistory();
                    } finally {
                        WaitCursorManager.removeWaitCursor(token);
                    }


                } else if (e.getKeyCode() == KeyEvent.VK_END) {
                    WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();
                    try {
                        frame.shiftOriginPixels(getWidth());
                        frame.recordHistory();
                    } finally {
                        WaitCursorManager.removeWaitCursor(token);
                    }
                } else if (e.getKeyCode() == KeyEvent.VK_PLUS) {
                } else if (e.getKeyCode() == KeyEvent.VK_MINUS) {
                }


            }
        };
        addKeyListener(keyAdapter);


        // Mouse Events
        MouseInputAdapter mouseAdapter = new DataPanelMouseAdapter();

        addMouseMotionListener(mouseAdapter);
        addMouseListener(mouseAdapter);
    }


    /**
     * Some performance hings from the Mac developer mailing list.  Might improve
     * graphics performance?
     * <p/>
     * // TODO  do timing tests with and without these hints
     *
     * @param g2D
     */
    private void applyMacPerformanceHints(Graphics2D g2D) {

        // From mac mailint list.  Might help performance ?
        g2D.setRenderingHint(RenderingHints.KEY_ALPHA_INTERPOLATION, RenderingHints.VALUE_ALPHA_INTERPOLATION_SPEED);
        g2D.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF);
        g2D.setRenderingHint(RenderingHints.KEY_COLOR_RENDERING, RenderingHints.VALUE_COLOR_RENDER_SPEED);
        g2D.setRenderingHint(RenderingHints.KEY_DITHERING, RenderingHints.VALUE_DITHER_DISABLE);
        g2D.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_NEAREST_NEIGHBOR);
        g2D.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_SPEED);


    }

    protected void addMousableRegion(MouseableRegion region) {
        parent.addMousableRegion(region);
    }

    protected void removeMousableRegions() {
        parent.getTrackRegions().clear();
    }

    public ReferenceFrame getFrame() {
        return frame;
    }


    /**
     * Receives all mouse events for a data panel.  Handling of some events are delegated to the current tool or track.
     */
    class DataPanelMouseAdapter extends MouseInputAdapter {

        /**
         * A scheduler is used to distinguish a click from a double click.
         */
        private ClickTaskScheduler clickScheduler = new ClickTaskScheduler();


        @Override
        public void mouseEntered(MouseEvent e) {
        }

        @Override
        public void mouseExited(MouseEvent e) {
            IGVMainFrame.getInstance().setSelectedRegion(null);
        }


        @Override
        public void mouseMoved(MouseEvent e) {
            tooltipTextPosition = e.getPoint();
            if (!frame.getChrName().equals(Globals.CHR_ALL)) {
                int location = (int) frame.getChromosomePosition(e.getX()) + 1;
                String position = frame.getChrName() + ":" + locationFormatter.format(location);
                IGVMainFrame.getInstance().setStatusBarMessage(position);
            }
        }

        /**
         * The mouse has been pressed.  If this is the platform's popup trigger select the track and popup a menu.
         * Otherwise delegate handling to the current tool.
         */
        @Override
        public void mousePressed(final MouseEvent e) {

            if (SwingUtilities.getWindowAncestor(DataPanel.this).isActive()) {
                DataPanel.this.requestFocus();
            }

            if (e.isPopupTrigger()) {
                IGVMainFrame.getInstance().getTrackManager().clearSelections();
                parent.selectTracks(e);
                TrackClickEvent te = new TrackClickEvent(e, frame);
                parent.openPopupMenu(te);
            } else {
                if (currentTool != null)
                    currentTool.mousePressed(e);
            }

        }

        /**
         * The mouse has been released.  If this is the platform's popup trigger select the track and popup a menu.
         * Otherwise delegate handling to the current tool.
         */
        @Override
        public void mouseReleased(MouseEvent e) {

            if (e.isPopupTrigger()) {
                IGVMainFrame.getInstance().getTrackManager().clearSelections();
                parent.selectTracks(e);
                TrackClickEvent te = new TrackClickEvent(e, frame);
                parent.openPopupMenu(te);
            } else {
                if (currentTool != null)
                    currentTool.mouseReleased(e);
            }
        }

        @Override
        public void mouseDragged(MouseEvent e) {
            if (currentTool != null)
                currentTool.mouseDragged(e);
        }

        @Override
        public void mouseClicked(final MouseEvent e) {
            if (e.getButton() != MouseEvent.BUTTON1 || e.isPopupTrigger()) {
                return;
            }

            // If this is the second click of a double click, cancel the scheduled single click task.
            // The shift and alt keys are alternative  zoom options
            //      shift zooms in by 8x,  alt zooms out by 2x
            if (e.getClickCount() > 1 || e.isShiftDown() || e.isAltDown() || (e.getClickCount() > 1)) {
                clickScheduler.cancelClickTask();

                int currentZoom = frame.getZoom();
                final int newZoom = e.isAltDown()
                        ? Math.max(currentZoom - 1, 0)
                        : (e.isShiftDown() ? currentZoom + 3 : currentZoom + 1);
                final double locationClicked = frame.getChromosomePosition(e.getX());

                WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();
                try {
                    frame.zoomTo(newZoom, locationClicked);
                } finally {
                    WaitCursorManager.removeWaitCursor(token);
                }
            } else if (e.getButton() == MouseEvent.BUTTON1 && !e.isPopupTrigger() && !e.isControlDown()) {

                // Unhandled single click.  Delegate to track unless second click arrives within double-click interval.
                // If the track does not handle the event delegate to the current tool
                TimerTask clickTask = new TimerTask() {

                    @Override
                    public void run() {
                        Object source = e.getSource();
                        if (source instanceof DataPanel) {
                            Track track = ((DataPanel) source).getTrack(e.getX(), e.getY());
                            if (track != null) {
                                TrackClickEvent te = new TrackClickEvent(e, frame);
                                if (track.handleDataClick(te)) {
                                    return;
                                } else {
                                    if (currentTool != null)
                                        currentTool.mouseClicked(e);
                                }
                            }
                        }
                    }

                };
                clickScheduler.scheduleClickTask(clickTask);
            }

        }
    }

}
