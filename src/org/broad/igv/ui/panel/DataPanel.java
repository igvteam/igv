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
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackClickEvent;
import org.broad.igv.track.TrackGroup;
import org.broad.igv.ui.*;
import org.broad.igv.ui.util.UIUtilities;

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

    private PanAndZoomTool panAndZoomTool;
    private IGVTool currentTool;
    private IGVTool defaultTool;
    private Point tooltipTextPosition;

    private ReferenceFrame frame;
    DataPanelContainer parent;
    private DataPanelPainter painter;


    public DataPanel(ReferenceFrame frame, DataPanelContainer parent) {
        init();
        this.frame = frame;
        this.parent = parent;
        this.setBackground(Color.white);
        setFocusable(true);
        setAutoscrolls(true);
        setToolTipText("Data panel");
        painter = new DataPanelPainter();
    }


    /**
     * @return the panAndZoomTool
     */
    public PanAndZoomTool getPanAndZoomTool() {
        return panAndZoomTool;
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


    /**
     * Method description
     *
     * @param tool
     */
    public void setCurrentTool(final IGVTool tool) {

        UIUtilities.invokeOnEventThread(new Runnable() {

            public void run() {

                if (defaultTool == null) {
                    defaultTool = tool;
                }

                // Swap out old tool
                if (currentTool != null) {
                    removeKeyListener(currentTool);
                    removeMouseListener(currentTool);
                    removeMouseMotionListener(currentTool);
                }

                // Swap in new tool
                currentTool = ((tool == null) ? defaultTool : tool);
                if (currentTool != null) {
                    setCursor(currentTool.getCursor());
                    addKeyListener(currentTool);
                    addMouseListener(currentTool);
                    addMouseMotionListener(currentTool);
                }
            }
        });
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

            context = new RenderContext(genomeId, this, graphics2D, frame);

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

            context = new RenderContext(genomeId, null, g, frame);
            final Collection<TrackGroup> groups = new ArrayList(parent.getTrackGroups());
            int trackWidth = rect.width;
            int trackHeight = rect.height;
            painter.paint(groups, context, trackWidth, trackHeight, getBackground(), rect);

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
        StringBuffer popupTextBuffer = new StringBuffer();
        double location = frame.getChromosomePosition(x);
        double displayLocation = location + 1;

        popupTextBuffer.append("<html>");

        Track track = null;
        List<MouseableRegion> regions = parent.getTrackRegions();
        for (MouseableRegion mouseRegion : regions) {
            if (mouseRegion.containsPoint(x, y)) {
                track = mouseRegion.getTracks().iterator().next();
                if (track != null) {

                    // Calculate y relative to track bounds
                    int yRelative = y - mouseRegion.getBounds().y;

                    // First see if there is an overlay track.  If there is, give
                    // it first crack
                    List<Track> overlays = IGVMainFrame.getInstance().getTrackManager().getOverlayTracks(track);

                    if (overlays != null) {
                        for (Track overlay : overlays) {
                            if ((overlay != track) && (overlay.getValueStringAt(
                                    frame.getChrName(), displayLocation, yRelative, frame) != null)) {
                                popupTextBuffer.append(getPopUpText(overlay, displayLocation, yRelative));
                                popupTextBuffer.append("<br>");
                            }
                        }
                    }

                    popupTextBuffer.append(getPopUpText(track, displayLocation, y));
                }
                break;
            }
        }

        popupTextBuffer.append("<br>Location: ");
        popupTextBuffer.append(frame.getChrName() + ":" + locationFormatter.format((int) displayLocation));

        setToolTipText(popupTextBuffer.toString());
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
     * Return the
     *
     * @return descriptive text for the panel at the given location.
     *         Note:  the coordinates are in panel coordinates, they must be
     *         offset by the track retangle.  This will not be needed when each
     *         track has its own panel.
     */
    String getPopUpText(Track track, double location, int y) {

        DataRange axisDefinition = track.getDataRange();
        StringBuffer buf = new StringBuffer();
        String value = track.getValueStringAt(frame.getChrName(), location, y, frame);
        if (value != null) {
            buf.append(value);
        }
        return buf.toString();
    }


    // TODO -- this method and the key and mouse adapters are nearly identical
    // for the DataPanel and FeaturePanel classes.  Refractor to combine those
    // classes, or share this method in some other way.

    private void init() {
        panAndZoomTool = new PanAndZoomTool(this);
        setCurrentTool(getPanAndZoomTool());
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


    class DataPanelMouseAdapter extends MouseInputAdapter {

        // TODO -- a static member on the frame class is probably not
        // the best place for this.

        @Override
        public void mouseEntered(MouseEvent e) {
        }

        @Override
        public void mouseExited(MouseEvent e) {
            IGVMainFrame.getInstance().setSelectedRegion(null);
        }

        /**
         * The mouse has been clicked.  If the mode is ZOOM_IN or ZOOM_OUT
         * zoom and center on click.  Otherwise ignore
         */
        @Override
        public void mousePressed(final MouseEvent e) {

            if (SwingUtilities.getWindowAncestor(DataPanel.this).isActive()) {
                DataPanel.this.requestFocus();
                //DataPanel.this.grabFocus();
            }

            if (e.isPopupTrigger()) {
                IGVMainFrame.getInstance().getTrackManager().clearSelections();
                parent.selectTracks(e);

                TrackClickEvent te = new TrackClickEvent(e, frame);
                parent.openPopupMenu(te);
            }
        }

        @Override
        public void mouseMoved(MouseEvent e) {
            tooltipTextPosition = e.getPoint();

        }

        @Override
        public void mouseReleased(MouseEvent e) {

            if (e.isPopupTrigger()) {
                IGVMainFrame.getInstance().getTrackManager().clearSelections();
                parent.selectTracks(e);
                TrackClickEvent te = new TrackClickEvent(e, frame);
                parent.openPopupMenu(te);
            }

        }
    }

}
