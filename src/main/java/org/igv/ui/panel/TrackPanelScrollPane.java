package org.igv.ui.panel;

import org.igv.Globals;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.ui.util.SnapshotUtilities;

import javax.swing.*;
import java.awt.*;

/**
 * A panel that wraps a track in a scrollable view with a fixed left strip
 * containing the drag handle and selection checkbox. The left strip stays
 * visible regardless of vertical scroll position.
 * <p>
 * Layout (null layout, manual positioning):
 * <pre>
 *   [SelectionPanel][DragHandle][ innerScrollPane (viewport → TrackPanel) ]
 * </pre>
 */
public class TrackPanelScrollPane extends JPanel implements Paintable {

    private static Logger log = LogManager.getLogger(TrackPanelScrollPane.class);

    TrackPanel trackPanel;
    private JScrollPane innerScrollPane;
    private DragHandlePanel dragHandlePanel;
    private TrackSelectionPanel selectionPanel;

    public TrackPanelScrollPane() {
        setLayout(null);
        setBorder(BorderFactory.createEmptyBorder());

        // Create the inner scroll pane for the track content
        innerScrollPane = new JScrollPane();
        innerScrollPane.setBorder(BorderFactory.createEmptyBorder());
        innerScrollPane.setForeground(new Color(153, 153, 153));
        innerScrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
        innerScrollPane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_NEVER);
        innerScrollPane.getVerticalScrollBar().setUnitIncrement(16);

        // Repaint name panels on end of scroll to center the text in the visible window.
        innerScrollPane.getVerticalScrollBar().addAdjustmentListener(e -> {
            if (!e.getValueIsAdjusting() && trackPanel != null) {
                trackPanel.getNamePanel().repaint();
            }
        });

        add(innerScrollPane);
    }

    // ---- Divider lines ----

    @Override
    protected void paintChildren(Graphics g) {
        super.paintChildren(g);
        drawPanelDividers(g);
    }

    /**
     * Draw vertical divider lines at the boundaries between name, attribute, and data panels.
     * Drawing at this level ensures divider lines survive when this scroll pane repaints
     * independently (e.g. after a track height change) without the parent container repainting.
     */
    private void drawPanelDividers(Graphics g) {
        if (trackPanel == null || trackPanel.mainPanel == null) return;
        MainPanel mainPanel = trackPanel.mainPanel;
        int leftOffset = mainPanel.getLeftOffset();
        int nameRight = leftOffset + mainPanel.getNamePanelX() + mainPanel.getNamePanelWidth();
        int dataLeft = leftOffset + mainPanel.getDataPanelX();

        int h = getHeight();
        Color dividerColor = Globals.isDarkMode() ? Color.GRAY : Color.LIGHT_GRAY;
        g.setColor(dividerColor);

        g.drawLine(nameRight, 0, nameRight, h);
        if (dataLeft != nameRight) {
            g.drawLine(dataLeft, 0, dataLeft, h);
        }
    }

    // ---- Layout ----

    @Override
    public void doLayout() {
        int h = getHeight();
        int w = getWidth();
        int selW = (selectionPanel != null && selectionPanel.isVisible())
                ? TrackSelectionPanel.SELECTION_PANEL_WIDTH : 0;
        int dragW = DragHandlePanel.DRAG_HANDLE_WIDTH;
        int leftW = selW + dragW;

        if (selectionPanel != null) {
            selectionPanel.setBounds(0, 0, selW, h);
        }
        if (dragHandlePanel != null) {
            dragHandlePanel.setBounds(selW, 0, dragW, h);
        }
        innerScrollPane.setBounds(leftW, 0, w - leftW, h);
    }

    // ---- Scroll pane delegation ----

    public JScrollBar getVerticalScrollBar() {
        return innerScrollPane.getVerticalScrollBar();
    }

    public int getVerticalScrollBarPolicy() {
        return innerScrollPane.getVerticalScrollBarPolicy();
    }

    public void setVerticalScrollBarPolicy(int policy) {
        innerScrollPane.setVerticalScrollBarPolicy(policy);
    }

    public JViewport getViewport() {
        return innerScrollPane.getViewport();
    }

    public Rectangle getViewportBorderBounds() {
        return innerScrollPane.getViewportBorderBounds();
    }

    // ---- Track panel management ----

    public void setViewportView(Component trackSetView) {
        if (!(trackSetView instanceof TrackPanel)) {
            throw new IllegalArgumentException("Class TrackPanelScrollPane can only contain a TrackPanel");
        }
        this.trackPanel = (TrackPanel) trackSetView;
        innerScrollPane.setViewportView(trackSetView);

        // Create left-strip components (fixed, outside the scrollable viewport)
        if (selectionPanel != null) remove(selectionPanel);
        if (dragHandlePanel != null) remove(dragHandlePanel);

        selectionPanel = new TrackSelectionPanel(trackPanel);
        selectionPanel.setVisible(TrackSelectionPanel.isSelectionModeActive());
        dragHandlePanel = new DragHandlePanel(trackPanel);

        add(selectionPanel);
        add(dragHandlePanel);

        updateScrollbarPolicy();
    }

    // ---- Sizing ----

    /**
     * Validates track height and triggers a full layout pass.
     */
    public void validateTrackHeight() {
        updateScrollbarPolicy();
        invalidate();
        Container parent = getParent();
        if (parent != null) {
            parent.invalidate();
            Component root = SwingUtilities.getRoot(this);
            if (root instanceof Container) {
                ((Container) root).validate();
            }
        }
        repaint();
    }

    /**
     * Updates the vertical scrollbar policy based on whether the content height
     * exceeds the scroll pane's preferred height.
     */
    public void updateScrollbarPolicy() {
        if (trackPanel != null && trackPanel.getTrack() != null) {
            int contentHeight = trackPanel.getTrack().getContentHeight();
            int scrollPaneHeight = getPreferredSize().height;
            boolean needsScrolling = contentHeight > scrollPaneHeight;

            int newPolicy = needsScrolling
                    ? ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS
                    : ScrollPaneConstants.VERTICAL_SCROLLBAR_NEVER;

            if (innerScrollPane.getVerticalScrollBarPolicy() != newPolicy) {
                innerScrollPane.setVerticalScrollBarPolicy(newPolicy);
            }
        }
    }

    @Override
    public Dimension getPreferredSize() {
        if (trackPanel == null || trackPanel.getTrack() == null) {
            return new Dimension(0, 0);
        }
        int height = PreferencesManager.getPreferences().getAsBoolean(Constants.SHOW_SINGLE_TRACK_PANE_KEY) ?
                trackPanel.getTrack().getContentHeight() : trackPanel.getTrack().getHeight();
        return new Dimension(super.getPreferredSize().width, height);
    }

    @Override
    public Dimension getMaximumSize() {
        Dimension pref = getPreferredSize();
        return new Dimension(Integer.MAX_VALUE, pref.height);
    }

    // ---- Accessors ----

    public TrackPanel getTrackPanel() {
        return trackPanel;
    }

    public String getTrackPanelName() {
        return trackPanel.getName();
    }

    public DataPanelContainer getDataPanel() {
        return trackPanel.getDataPanelContainer();
    }

    public boolean isEmpty() {
        return !trackPanel.hasTracks();
    }

    public TrackNamePanel getNamePanel() {
        return trackPanel.getNamePanel();
    }

    public DragHandlePanel getDragHandlePanel() {
        return dragHandlePanel;
    }

    public TrackSelectionPanel getSelectionPanel() {
        return selectionPanel;
    }

    /**
     * Show or hide the selection panel (checkbox panel).
     */
    public void setSelectionPanelVisible(boolean visible) {
        TrackSelectionPanel.setSelectionModeActive(visible);
        if (selectionPanel != null) {
            if (!visible) {
                selectionPanel.setTrackSelected(false);
            }
            selectionPanel.setVisible(visible);
        }
        doLayout();
        revalidate();
        repaint();
    }

    /** Returns the width of the fixed left strip (drag handle + optional selection). */
    int getLeftStripWidth() {
        int w = DragHandlePanel.DRAG_HANDLE_WIDTH;
        if (selectionPanel != null && selectionPanel.isVisible()) {
            w += TrackSelectionPanel.SELECTION_PANEL_WIDTH;
        }
        return w;
    }

    // ---- Appearance ----

    @Override
    public void setBackground(Color color) {
        super.setBackground(color);
        if (innerScrollPane != null) innerScrollPane.setBackground(color);
        if (dragHandlePanel != null) dragHandlePanel.setBackground(color);
        if (selectionPanel != null) selectionPanel.setBackground(color);
        if (trackPanel != null) trackPanel.setBackground(color);
    }

    // ---- Offscreen rendering ----

    @Override
    public void paintOffscreen(Graphics2D g, Rectangle tspRect, boolean batch) {
        Graphics2D panelGraphics = (Graphics2D) g.create();
        // Shift right to account for the fixed left strip
        int leftW = getLeftStripWidth();
        panelGraphics.translate(leftW, 0);
        int yOffset = getViewport().getViewPosition().y;
        panelGraphics.translate(0, -yOffset);
        Rectangle panelRect = getViewportBorderBounds();
        panelRect.y = yOffset;
        panelRect.height = tspRect.height;
        trackPanel.paintOffscreen(panelGraphics, panelRect, batch);
        panelGraphics.dispose();
    }

    @Override
    public int getSnapshotHeight(boolean batch) {
        if (trackPanel.getTracks().isEmpty()) {
            return 0;
        }
        if (batch) {
            int panelHeight;
            int maxPanelHeight = SnapshotUtilities.getMaxPanelHeight();
            final int scrollPaneHeight = getHeight();
            if (maxPanelHeight <= 0) {
                panelHeight = scrollPaneHeight;
            } else {
                int contentHeight = trackPanel.getContentHeight();
                panelHeight = Math.min(maxPanelHeight, Math.max(scrollPaneHeight, contentHeight));
            }
            return panelHeight;
        } else {
            return getHeight();
        }
    }
}
