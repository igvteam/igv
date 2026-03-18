package org.igv.ui.panel;

import org.igv.Globals;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.track.Track;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

/**
 * A thin draggable divider between TrackPanelScrollPanes. Dragging the divider
 * resizes the panel above (expanding or shrinking) and the panel below
 * (shrinking or expanding correspondingly). Height changes are persisted via
 * {@link Track#setHeight(int)}.
 */
public class TrackPanelDivider extends JPanel {

    public static final int DIVIDER_HEIGHT = 5;
    private static final int MIN_TRACK_HEIGHT = 20;

    private final TrackPanelScrollPane abovePane;
    private final TrackPanelScrollPane belowPane;

    private int dragStartY;
    private int originalAboveHeight;
    private int originalBelowHeight;

    public TrackPanelDivider(TrackPanelScrollPane abovePane, TrackPanelScrollPane belowPane) {
        this.abovePane = abovePane;
        this.belowPane = belowPane;

        setCursor(Cursor.getPredefinedCursor(Cursor.N_RESIZE_CURSOR));
        setPreferredSize(new Dimension(Integer.MAX_VALUE, DIVIDER_HEIGHT));
        setMinimumSize(new Dimension(0, DIVIDER_HEIGHT));
        setMaximumSize(new Dimension(Integer.MAX_VALUE, DIVIDER_HEIGHT));
        updateBackground();

        MouseAdapter mouseAdapter = new MouseAdapter() {
            @Override
            public void mousePressed(MouseEvent e) {
                // Don't allow resizing in single-track-pane mode (content-driven heights)
                if (PreferencesManager.getPreferences().getAsBoolean(Constants.SHOW_SINGLE_TRACK_PANE_KEY)) {
                    return;
                }
                dragStartY = e.getYOnScreen();
                originalAboveHeight = getTrackHeight(abovePane);
                originalBelowHeight = getTrackHeight(belowPane);
            }

            @Override
            public void mouseDragged(MouseEvent e) {
                if (PreferencesManager.getPreferences().getAsBoolean(Constants.SHOW_SINGLE_TRACK_PANE_KEY)) {
                    return;
                }
                int delta = e.getYOnScreen() - dragStartY;

                int newAboveHeight = originalAboveHeight + delta;
                int newBelowHeight = originalBelowHeight - delta;

                // Clamp to minimum heights
                if (newAboveHeight < MIN_TRACK_HEIGHT) {
                    newBelowHeight += (MIN_TRACK_HEIGHT - newAboveHeight);
                    newAboveHeight = MIN_TRACK_HEIGHT;
                }
                if (newBelowHeight < MIN_TRACK_HEIGHT) {
                    newAboveHeight += (MIN_TRACK_HEIGHT - newBelowHeight);
                    newBelowHeight = MIN_TRACK_HEIGHT;
                }

                // Set the heights on the Track objects — this persists in sessions
                setTrackHeight(abovePane, newAboveHeight);
                setTrackHeight(belowPane, newBelowHeight);
            }
        };

        addMouseListener(mouseAdapter);
        addMouseMotionListener(mouseAdapter);
    }

    private static int getTrackHeight(TrackPanelScrollPane pane) {
        Track track = pane.getTrackPanel().getTrack();
        return track != null ? track.getHeight() : 0;
    }

    private static void setTrackHeight(TrackPanelScrollPane pane, int height) {
        Track track = pane.getTrackPanel().getTrack();
        if (track != null) {
            track.setHeight(height);
        }
    }

    private void updateBackground() {
        if (Globals.isDarkMode()) {
            setBackground(UIManager.getColor("Panel.background"));
        } else {
            setBackground(new Color(230, 230, 230));
        }
    }

    @Override
    protected void paintComponent(Graphics g) {
        updateBackground();
        super.paintComponent(g);

        Graphics2D g2d = (Graphics2D) g.create();
        int w = getWidth();
        int h = getHeight();
        int centerY = h / 2;

        // Draw a subtle grip line in the center
        g2d.setColor(Globals.isDarkMode() ? new Color(120, 120, 120) : new Color(180, 180, 180));
        int lineMargin = Math.max(4, w / 4);
        g2d.drawLine(lineMargin, centerY, w - lineMargin, centerY);

        g2d.dispose();
    }
}

