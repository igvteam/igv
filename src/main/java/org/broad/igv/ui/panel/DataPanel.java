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

import com.google.common.base.Objects;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.event.DataLoadedEvent;
import org.broad.igv.event.IGVEventObserver;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackClickEvent;
import org.broad.igv.track.TrackGroup;
import org.broad.igv.ui.AbstractDataPanelTool;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.ui.util.DataPanelTool;
import org.broad.igv.ui.util.MessageUtils;

import javax.swing.*;
import javax.swing.event.MouseInputAdapter;
import java.awt.*;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.text.DecimalFormat;
import java.util.List;
import java.util.*;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.Collectors;

/**
 * The batch panel for displaying tracks and data.  A DataPanel is always associated with a ReferenceFrame.  Normally
 * there is a single reference frame (and thus panel), but when "gene list" or other split screen views are
 * invoked there can be multiple panels.
 *
 * @author jrobinso
 */
public class DataPanel extends JComponent implements Paintable, IGVEventObserver {

    private static Logger log = Logger.getLogger(DataPanel.class);

    private DataPanelTool defaultTool;
    private DataPanelTool currentTool;
    private ReferenceFrame frame;
    private DataPanelContainer parent;
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
        setBackground(PreferencesManager.getPreferences().getAsColor(Constants.BACKGROUND_COLOR));
        ToolTipManager.sharedInstance().registerComponent(this);
    }

    @Override
    public void receiveEvent(Object event) {
        if (event instanceof DataLoadedEvent) {
            if (((DataLoadedEvent) event).referenceFrame == frame) {
                log.debug("Data loaded repaint " + frame);
                repaint();
            }
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

    long lastPaintTime = 0;

    @Override
    public void paintComponent(final Graphics g) {

        super.paintComponent(g);
        RenderContext context = null;
        try {

            lastPaintTime = System.currentTimeMillis();

            Rectangle clipBounds = g.getClipBounds();
            final Rectangle visibleRect = getVisibleRect();
            final Rectangle damageRect = clipBounds == null ? visibleRect : clipBounds.intersection(visibleRect);
            Graphics2D graphics2D = (Graphics2D) g; //(Graphics2D) g.create();

            context = new RenderContext(this, graphics2D, frame, visibleRect);

            final Collection<TrackGroup> groups = parent.getTrackGroups();

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


            long dt = System.currentTimeMillis() - lastPaintTime;
            PanTool.repaintTime(dt);

        } catch (Exception e) {
            MessageUtils.showMessage("Unexpected error repainting view.  See igv.log for details.");
            log.error("Error painting DataPanel", e);
        } finally {
            if (context != null) {
                context.dispose();
            }
        }
    }

    @Override
    public void setBounds(int x, int y, int width, int height) {
        super.setBounds(x, y, width, height);
        Insets insets = this.getInsets();
        frame.setBounds(x + insets.left, width - insets.left - insets.right);
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

                List<Track> trackList = group.getVisibleTracks();
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

    public void paintOffscreen(final Graphics2D g, Rectangle rect, boolean batch) {

        RenderContext context = null;
        Graphics borderGraphics = g.create();
        borderGraphics.setColor(Color.darkGray);

        try {

            context = new RenderContext(null, g, frame, rect);
            final Collection<TrackGroup> groups = new ArrayList(parent.getTrackGroups());
            int width = rect.width;
            painter.paint(groups, context, width, getBackground(), rect);

            drawAllRegions(g);

            borderGraphics.drawRect(0, rect.y, rect.width - 1, rect.height - 1);

        } finally {
            if (context != null) {
                context.dispose();
            }
            borderGraphics.dispose();
        }
    }

    @Override
    public int getSnapshotHeight(boolean batch) {
        return getHeight();
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

        boolean drawBars = PreferencesManager.getPreferences().getAsBoolean(Constants.SHOW_REGION_BARS);
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
                                    frame.getChrName(), position, x, y, frame) != null)) {
                                String valueString = overlay.getValueStringAt(frame.getChrName(), position, x, y, frame);
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
                        String valueString = track.getValueStringAt(frame.getChrName(), position, x, y, frame);
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
        setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(200, 200, 200)));
        setRequestFocusEnabled(false);

        // Key Events
        KeyAdapter keyAdapter = new KeyAdapter() {

            @Override
            public void keyPressed(KeyEvent e) {
                int shiftOriginPixels = Integer.MIN_VALUE;
                int zoomIncr = Integer.MIN_VALUE;
                boolean showWaitCursor = false;

                if (e.getKeyChar() == '+' || e.getKeyCode() == KeyEvent.VK_PLUS) {
                    zoomIncr = +1;
                    showWaitCursor = true;
                } else if (e.getKeyChar() == '-' || e.getKeyCode() == KeyEvent.VK_PLUS) {
                    zoomIncr = -1;
                    showWaitCursor = true;
                } else if (e.getKeyCode() == KeyEvent.VK_RIGHT) {
                    shiftOriginPixels = 50;
                } else if (e.getKeyCode() == KeyEvent.VK_LEFT) {
                    shiftOriginPixels = -50;
                } else if (e.getKeyCode() == KeyEvent.VK_HOME) {
                    shiftOriginPixels = -getWidth();
                    showWaitCursor = true;
                } else if (e.getKeyCode() == KeyEvent.VK_END) {
                    shiftOriginPixels = getWidth();
                    showWaitCursor = true;
                } else if (e.getKeyCode() == KeyEvent.VK_PLUS) {
                } else if (e.getKeyCode() == KeyEvent.VK_MINUS) {
                }

                WaitCursorManager.CursorToken token = null;
                if (showWaitCursor) token = WaitCursorManager.showWaitCursor();
                try {
                    if (zoomIncr > Integer.MIN_VALUE) {
                        frame.doZoomIncrement(zoomIncr);
                    } else if (shiftOriginPixels > Integer.MIN_VALUE) {
                        frame.shiftOriginPixels(shiftOriginPixels);
                    } else {
                        return;
                    }

                    //Assume that anything special enough to warrant a wait cursor
                    //should be in history
                    if (showWaitCursor) {
                        frame.recordHistory();
                    }
                } finally {
                    if (token != null) WaitCursorManager.removeWaitCursor(token);
                }

            }
        };
        addKeyListener(keyAdapter);


        // Mouse Events
        MouseInputAdapter mouseAdapter = new DataPanelMouseAdapter();

        addMouseMotionListener(mouseAdapter);
        addMouseListener(mouseAdapter);
        addMouseWheelListener(mouseAdapter);
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

        long lastClickTime = 0;
        MouseEvent mouseDown = null;


        @Override
        public void mouseMoved(MouseEvent e) {
            String position = null;
            if (!frame.getChrName().equals(Globals.CHR_ALL)) {
                int location = (int) frame.getChromosomePosition(e.getX()) + 1;
                position = frame.getChrName() + ":" + locationFormatter.format(location);
                IGV.getInstance().setStatusBarMessag2(position);
            }
            updateTooltipText(e.getX(), e.getY());

            if (IGV.getInstance().isRulerEnabled()) {
                IGV.getInstance().repaint();
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
                doPopupMenu(e);
                e.consume();
            } else {
                if (currentTool != null)
                    currentTool.mousePressed(e);
                mouseDown = e;
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
                e.consume();
            } else {
                if (mouseDown != null && distance(mouseDown, e) < 5) {
                    doMouseClick(e);
                } else if (currentTool != null)
                    currentTool.mouseReleased(e);
            }
            mouseDown = null;
        }

        private void doPopupMenu(MouseEvent e) {
            IGV.getInstance().clearSelections();
            parent.selectTracks(e);
            TrackClickEvent te = new TrackClickEvent(e, frame);
            parent.openPopupMenu(te);
        }

        private double distance(MouseEvent e1, MouseEvent e2) {
            double dx = e1.getX() - e2.getX();
            double dy = e1.getY() - e2.getY();
            return Math.sqrt(dx * dx + dy * dy);
        }


        /**
         * The mouse has been dragged.  Delegate to current tool.
         *
         * @param e
         */
        @Override
        public void mouseDragged(MouseEvent e) {
            if (mouseDown != null && currentTool != null)
                currentTool.mouseDragged(e);
        }

        @Override
        public void mouseEntered(MouseEvent e) {
            mouseDown = null;
        }

        @Override
        public void mouseExited(MouseEvent e) {
            mouseDown = null;
        }

        /**
         * Zoom in/out when modifier + scroll wheel used
         *
         * @param e
         */
        @Override
        public void mouseWheelMoved(MouseWheelEvent e) {
            //we use either ctrl or meta to deal with PCs and Macs
            if (e.isControlDown() || e.isMetaDown()) {
                int wheelRotation = e.getWheelRotation();
                //Mouse move up is negative, that should zoom in
                int zoomIncr = -wheelRotation / 2;
                getFrame().doZoomIncrement(zoomIncr);
            }
            //TODO Use this to pan. Seems weird, but it's how side scrolling on my mouse gets interpreted,
            //so could be handy for people with 2D wheels
//            else if(e.isShiftDown()){
//                System.out.println(e);
//            }
            else {
                //Default action if no modifier
                e.getComponent().getParent().dispatchEvent(e);
            }
        }


        /**
         * The mouse was clicked. If this is the second click of a double click, cancel the scheduled single click task.
         * The shift and alt keys are alternative  zoom options
         * shift zooms in by 8x,  alt zooms out by 2x
         * <p>
         * NOTE: mouseClick is not used because in Java a mouseClick event is emitted only if the mouse has not
         * moved at all between press and release.  This is difficult to do, even when trying.
         * <p>
         * <p/>
         * TODO -- the "currentTool" is also a mouselistener, so there are two.  This makes mouse event handling
         * TODO -- needlessly complicated, which handler has preference, etc.  Move this code to the default
         * TODO -- PanAndZoomTool
         *
         * @param e
         */

        public void doMouseClick(final MouseEvent e) {

            long clickTime = System.currentTimeMillis();

            if (e.isPopupTrigger()) {
                return;
            }

            if (currentTool instanceof RegionOfInterestTool) {
                currentTool.mouseClicked(e);
                e.consume();
                return;
            }

            if (e.isPopupTrigger()) {
                doPopupMenu(e);
                e.consume();
                return;
            }

            Object source = e.getSource();
            if (source instanceof DataPanel && e.getButton() == MouseEvent.BUTTON1) {
                final Track track = ((DataPanel) e.getSource()).getTrack(e.getX(), e.getY());

                if (e.isShiftDown()) {
                    final double locationClicked = frame.getChromosomePosition(e.getX());
                    frame.doIncrementZoom(3, locationClicked);
                    e.consume();
                } else if (e.isAltDown()) {
                    final double locationClicked = frame.getChromosomePosition(e.getX());
                    frame.doIncrementZoom(-1, locationClicked);
                    e.consume();
                } else if ((e.isMetaDown() || e.isControlDown()) && track != null) {
                    TrackClickEvent te = new TrackClickEvent(e, frame);
                    if (track.handleDataClick(te)) {
                        e.consume();
                        return;
                    }

                } else {

                    // No modifier, left-click.  Defer processing with a timer until we are sure this is not the
                    // first of a "double-click".

                    if (clickTime - lastClickTime < UIConstants.getDoubleClickInterval()) {
                        clickScheduler.cancelClickTask();
                        final double locationClicked = frame.getChromosomePosition(e.getX());
                        frame.doIncrementZoom(1, locationClicked);

                    } else {

                        lastClickTime = clickTime;

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
