/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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

import com.google.common.base.Objects;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.track.*;
import org.broad.igv.ui.AbstractDataPanelTool;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.ui.util.DataPanelTool;

import javax.swing.*;
import javax.swing.event.MouseInputAdapter;
import java.awt.*;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.text.DecimalFormat;
import java.util.*;
import java.util.List;

/**
 * The batch panel for displaying tracks and data.  A DataPanel is always associated with a ReferenceFrame.  Normally
 * there is a single reference frame (and thus panel), but when "gene list" or other split screen views are
 * invoked there can be multiple panels.
 *
 * @author jrobinso
 */
public class DataPanel extends JComponent implements Paintable {

    private static Logger log = Logger.getLogger(DataPanel.class);
    private boolean isWaitingForToolTipText = false;

    private DataPanelTool defaultTool;
    private DataPanelTool currentTool;
    // private Point tooltipTextPosition;
    private ReferenceFrame frame;
    DataPanelContainer parent;
    private DataPanelPainter painter;
    private String tooltipText = "";


    public DataPanel(ReferenceFrame frame, DataPanelContainer parent) {
        init();
        this.defaultTool = new PanTool(this);
        this.currentTool = defaultTool;
        this.frame = frame;
        this.parent = parent;
        setFocusable(true);
        setAutoscrolls(true);
        setToolTipText("");
        painter = new DataPanelPainter();
        setBackground(PreferenceManager.getInstance().getAsColor(PreferenceManager.BACKGROUND_COLOR));

        ToolTipManager.sharedInstance().registerComponent(this);
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
        ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        RenderContext context = null;
        try {

            Rectangle clipBounds = g.getClipBounds();
            final Rectangle visibleRect = getVisibleRect();
            final Rectangle damageRect = clipBounds == null ? visibleRect : clipBounds.intersection(visibleRect);
            Graphics2D graphics2D = (Graphics2D) g; //(Graphics2D) g.create();

            context = new RenderContextImpl(this, graphics2D, frame, visibleRect);

            if (Globals.IS_MAC) {
                this.applyMacPerformanceHints((Graphics2D) g);
            }

            final Collection<TrackGroup> groups = parent.getTrackGroups();

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

            computeMousableRegions(groups, trackWidth);
            painter.paint(groups, context, trackWidth, getBackground(), damageRect);


            // If there is a partial ROI in progress draw it first
            if (currentTool instanceof RegionOfInterestTool) {
                int startLoc = ((RegionOfInterestTool) currentTool).getRoiStart();
                if (startLoc > 0) {
                    int start = frame.getScreenPosition(startLoc);
                    g.setColor(Color.BLACK);
                    graphics2D.drawLine(start, 0, start, getHeight());
                }
            }

            drawAllRegions(g);


        } finally {

            if (context != null) {
                context.dispose();
            }
        }
    }

    /**
     * TODO -- move this to a "layout" command, to layout tracks and assign positions
     */
    private void computeMousableRegions(Collection<TrackGroup> groups, int width) {

        final List<MouseableRegion> mouseableRegions = parent.getMouseRegions();
        mouseableRegions.clear();
        int trackX = 0;
        int trackY = 0;
        for (Iterator<TrackGroup> groupIter = groups.iterator(); groupIter.hasNext(); ) {
            TrackGroup group = groupIter.next();


            if (group.isVisible()) {
                if (groups.size() > 1) {
                    trackY += UIConstants.groupGap;
                }

                List<Track> trackList = group.getTracks();
                for (Track track : trackList) {
                    if (track == null) continue;
                    int trackHeight = track.getHeight();

                    if (track.isVisible()) {
                        Rectangle rect = new Rectangle(trackX, trackY, width, trackHeight);
                        if (mouseableRegions != null) {
                            mouseableRegions.add(new MouseableRegion(rect, track));
                        }
                        trackY += trackHeight;
                    }
                }

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

            context = new RenderContextImpl(null, g, frame, rect);
            final Collection<TrackGroup> groups = new ArrayList(parent.getTrackGroups());
            int width = rect.width;
            painter.paint(groups, context, width, getBackground(), rect);

            drawAllRegions(g);

            Color c = g.getColor();
            g.setColor(Color.darkGray);
            g.drawRect(rect.x, rect.y, rect.width, rect.height);
            g.setColor(c);            //super.paintBorder(g);

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

        // TODO -- get rid of this ugly reference to IGV
        Collection<RegionOfInterest> regions =
                IGV.getInstance().getSession().getRegionsOfInterest(frame.getChrName());

        if ((regions == null) || regions.isEmpty()) {
            return;
        }

        boolean drawBars = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_REGION_BARS);
        Graphics2D graphics2D = (Graphics2D) g.create();
        graphics2D.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
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
        for (MouseableRegion mouseRegion : parent.getMouseRegions()) {
            if (mouseRegion.containsPoint(x, y)) {
                return mouseRegion.getTracks().iterator().next();
            }
        }
        return null;

    }


    @Override
    public void setToolTipText(String text) {
        if (!Objects.equal(tooltipText, text)) {
            IGV.getInstance().setStatusWindowText(text);
            this.tooltipText = text;
            putClientProperty(TOOL_TIP_TEXT_KEY, text);
        }

    }

    /**
     * {@inheritDoc}
     * <p/>
     * The tooltip text may be null, in which case no tooltip is displayed
     */
    @Override
    final public String getToolTipText() {
        //TODO Suppress tooltips instead. This is hard to get exactly right
        //TODO with our different tooltip settings
        if (currentTool instanceof RegionOfInterestTool) {
            return null;
        }
        return tooltipText;
    }


    /**
     * Update tooltip text for the current mouse position (x, y)
     *
     * @param x Mouse x position in pixels
     * @param y Mouse y position in pixels
     */
    public void updateTooltipText(int x, int y) {

        //Tooltip here specifically means text that is shown on hover
        //We disable it unless that option is specified
        if (!IGV.getInstance().isShowDetailsOnHover()) {
            setToolTipText(null);
            return;
        }

        double position = frame.getChromosomePosition(x);

        Track track = null;
        List<MouseableRegion> regions = parent.getMouseRegions();
        StringBuffer popupTextBuffer = new StringBuffer();
        popupTextBuffer.append("<html>");

        for (MouseableRegion mouseRegion : regions) {
            if (mouseRegion.containsPoint(x, y)) {
                track = mouseRegion.getTracks().iterator().next();
                if (track != null) {

                    // First see if there is an overlay track.  If there is, give
                    // it first crack
                    List<Track> overlays = IGV.getInstance().getOverlayTracks(track);
                    boolean foundOverlaidFeature = false;
                    if (overlays != null) {
                        for (Track overlay : overlays) {
                            if ((overlay != track) && (overlay.getValueStringAt(
                                    frame.getChrName(), position, y, frame) != null)) {
                                String valueString = overlay.getValueStringAt(frame.getChrName(), position, y, frame);
                                if (valueString != null) {
                                    popupTextBuffer.append(valueString);
                                    popupTextBuffer.append("<br>");
                                    foundOverlaidFeature = true;
                                    break;
                                }
                            }
                        }
                    }
                    if (!foundOverlaidFeature) {
                        String valueString = track.getValueStringAt(frame.getChrName(), position, y, frame);
                        if (valueString != null) {
                            if (foundOverlaidFeature) {
                                popupTextBuffer.append("---------------------<br>");
                            }
                            popupTextBuffer.append(valueString);
                            popupTextBuffer.append("<br>");
                            break;
                        }
                    }
                }
            }
        }

        if (popupTextBuffer.length() > 6) {   // 6 characters for <html>
            //popupTextBuffer.append("<br>--------------------------");
            //popupTextBuffer.append(positionString);
            String puText = popupTextBuffer.toString().trim();
            if (!puText.equals(tooltipText)) {
                setToolTipText(puText);
            }
        } else {
            setToolTipText(null);
        }
    }


    private void init() {
        setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
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

    protected void removeMousableRegions() {
        parent.getMouseRegions().clear();
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
        public void mouseMoved(MouseEvent e) {
            String position = null;
            if (!frame.getChrName().equals(Globals.CHR_ALL)) {
                int location = (int) frame.getChromosomePosition(e.getX()) + 1;
                position = frame.getChrName() + ":" + locationFormatter.format(location);
                IGV.getInstance().setStatusBarPosition(position);
            }
            updateTooltipText(e.getX(), e.getY());

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
                doPopupMenu(e);
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
                doPopupMenu(e);
            } else {
                if (currentTool != null)
                    currentTool.mouseReleased(e);
            }
        }

        private void doPopupMenu(MouseEvent e) {
            IGV.getInstance().clearSelections();
            parent.selectTracks(e);
            TrackClickEvent te = new TrackClickEvent(e, frame);
            parent.openPopupMenu(te);
        }


        /**
         * The mouse has been dragged.  Delegate to current tool.
         *
         * @param e
         */
        @Override
        public void mouseDragged(MouseEvent e) {
            if (currentTool != null)
                currentTool.mouseDragged(e);
        }


        /**
         * The mouse was clicked. If this is the second click of a double click, cancel the scheduled single click task.
         * The shift and alt keys are alternative  zoom options
         * shift zooms in by 8x,  alt zooms out by 2x
         * <p/>
         * TODO -- the "currentTool" is also a mouselistener, so there are two.  This makes mouse event handling
         * TODO -- needlessly complicated, which handler has preference, etc.  Move this code to the default
         * TODO -- PanAndZoomTool
         *
         * @param e
         */
        @Override
        public void mouseClicked(final MouseEvent e) {

            // ctrl-mouse down is the mac popup trigger, but you will also get a clck even.  Ignore the click.
            if (Globals.IS_MAC && e.isControlDown()) {
                return;
            }

            if (currentTool instanceof RegionOfInterestTool) {
                currentTool.mouseClicked(e);
                return;
            }

            if (e.isPopupTrigger()) {
                doPopupMenu(e);
                return;
            }

            Object source = e.getSource();
            if (source instanceof DataPanel && e.getButton() == MouseEvent.BUTTON1) {
                final Track track = ((DataPanel) e.getSource()).getTrack(e.getX(), e.getY());

                if (e.isShiftDown()) {
                    final double locationClicked = frame.getChromosomePosition(e.getX());
                    frame.zoomBy(3, locationClicked);
                } else if (e.isAltDown()) {
                    final double locationClicked = frame.getChromosomePosition(e.getX());
                    frame.zoomBy(-1, locationClicked);
                } else if ((e.isMetaDown() || e.isControlDown()) && track != null) {
                    TrackClickEvent te = new TrackClickEvent(e, frame);

                    track.handleDataClick(te);

                } else {

                    // No modifier, left-click.  Defer processing with a timer until we are sure this is not the
                    // first of a "double-click".

                    if (e.getClickCount() > 1) {
                        clickScheduler.cancelClickTask();
                        final double locationClicked = frame.getChromosomePosition(e.getX());
                        frame.zoomBy(1, locationClicked);

                    } else {

                        // Unhandled single click.  Delegate to track or tool unless second click arrives within
                        // double-click interval.
                        TimerTask clickTask = new TimerTask() {
                            @Override
                            public void run() {
                                Object source = e.getSource();
                                if (source instanceof DataPanel) {

                                    if (track != null) {
                                        TrackClickEvent te = new TrackClickEvent(e, frame);
                                        List<Track> overlays = IGV.getInstance().getOverlayTracks(track);
                                        boolean handled = false;
                                        if (overlays != null) {
                                            for (Track overlay : overlays) {
                                                if (overlay.getFeatureAtMousePosition(te) != null) {
                                                    overlay.handleDataClick(te);
                                                    handled = true;
                                                }
                                            }
                                        }
                                        if (!handled) {
                                            handled = track.handleDataClick(te);
                                        }


                                        if (handled) {
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
    }

    /**
     * A utility class for sceduling single-click actions "in the future",
     *
     * @author jrobinso
     * @date Dec 17, 2010
     */
    public class PopupTextUpdater {

        private TimerTask currentClickTask;

        public void cancelClickTask() {
            if (currentClickTask != null) {
                currentClickTask.cancel();
                currentClickTask = null;
            }
        }

        public void scheduleUpdateTask(TimerTask task) {
            cancelClickTask();
            currentClickTask = task;
            (new java.util.Timer()).schedule(currentClickTask, UIConstants.getDoubleClickInterval());
        }
    }
}