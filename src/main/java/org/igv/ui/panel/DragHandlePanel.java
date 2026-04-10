package org.igv.ui.panel;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

/**
 * A drag handle panel that allows users to reorder TrackPanels within their container
 * by dragging. Displays a grip pattern to indicate draggability.
 */
public class DragHandlePanel extends JPanel {

    public static final int DRAG_HANDLE_WIDTH = 12;

    private final TrackPanel trackPanel;
    private boolean isDragging = false;

    public DragHandlePanel(TrackPanel trackPanel) {
        this.trackPanel = trackPanel;
        setBackground(Color.WHITE);
        setCursor(Cursor.getPredefinedCursor(Cursor.MOVE_CURSOR));
        setPreferredSize(new Dimension(DRAG_HANDLE_WIDTH, 0));
        setMinimumSize(new Dimension(DRAG_HANDLE_WIDTH, 0));

        MouseAdapter mouseAdapter = new MouseAdapter() {
            @Override
            public void mousePressed(MouseEvent e) {
                isDragging = false;
            }

            @Override
            public void mouseReleased(MouseEvent e) {
                isDragging = false;
            }

            @Override
            public void mouseDragged(MouseEvent e) {
                if (isDragging) {
                    return;
                }
                isDragging = true;
                JComponent c = trackPanel;
                TransferHandler handler = c.getTransferHandler();
                if (handler != null) {
                    handler.exportAsDrag(c, e, TransferHandler.MOVE);
                }
            }
        };

        addMouseListener(mouseAdapter);
        addMouseMotionListener(mouseAdapter);
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);

        Graphics2D g2d = (Graphics2D) g.create();
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        // DragHandlePanel is outside the scrollable viewport, so getHeight()
        // directly gives the visible height — no getVisibleRect() workaround needed.
        int width = getWidth();
        int panelHeight = getHeight();

        g2d.setColor(new Color(180, 180, 180));

        // Draw grip pattern centered vertically
        int gripHeight = Math.min(24, panelHeight - 8);
        if (gripHeight > 0) {
            int startY = (panelHeight - gripHeight) / 2;
            int lineSpacing = 4;
            int lineWidth = width - 4;
            int startX = 2;

            for (int y = startY; y < startY + gripHeight; y += lineSpacing) {
                g2d.drawLine(startX, y, startX + lineWidth, y);
            }
        }

        g2d.dispose();
    }

    public TrackPanel getTrackPanel() {
        return trackPanel;
    }
}

