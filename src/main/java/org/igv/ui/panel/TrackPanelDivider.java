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
 * A thin draggable divider placed below a TrackPanelScrollPane. Dragging works
 * like Google Sheets row resizing:
 * <ul>
 *   <li><b>Drag down</b> – increases the height of the track above; everything
 *       below simply shifts down.</li>
 *   <li><b>Drag up</b> – decreases the height of the track above (down to
 *       {@link Track#getMinimumHeight()}); everything below simply shifts up.</li>
 * </ul>
 * Only the track above the divider is resized; tracks below are never modified.
 * Height changes are persisted via {@link Track#setHeight(int)}.
 */
public class TrackPanelDivider extends JPanel {

    public static final int DIVIDER_HEIGHT = 5;

    /** The pane above this divider (whose height shrinks/grows). */
    private final TrackPanelScrollPane abovePane;

    private int dragStartY;
    private int originalAboveHeight;

    /**
     * @param abovePane the scroll pane above this divider, or {@code null} if none
     * @param belowPane the scroll pane below this divider, or {@code null} for the trailing divider
     */
    public TrackPanelDivider(TrackPanelScrollPane abovePane, TrackPanelScrollPane belowPane) {
        this.abovePane = abovePane;

        setCursor(Cursor.getPredefinedCursor(Cursor.N_RESIZE_CURSOR));
        setPreferredSize(new Dimension(Integer.MAX_VALUE, DIVIDER_HEIGHT));
        setMinimumSize(new Dimension(0, DIVIDER_HEIGHT));
        setMaximumSize(new Dimension(Integer.MAX_VALUE, DIVIDER_HEIGHT));
        updateBackground();

        MouseAdapter mouseAdapter = new MouseAdapter() {
            @Override
            public void mousePressed(MouseEvent e) {
                if (abovePane == null) return;
                if (PreferencesManager.getPreferences().getAsBoolean(Constants.SHOW_SINGLE_TRACK_PANE_KEY)) {
                    return;
                }
                dragStartY = e.getYOnScreen();
                originalAboveHeight = getTrackHeight(abovePane);
            }

            @Override
            public void mouseDragged(MouseEvent e) {
                if (abovePane == null) return;
                if (PreferencesManager.getPreferences().getAsBoolean(Constants.SHOW_SINGLE_TRACK_PANE_KEY)) {
                    return;
                }
                int delta = e.getYOnScreen() - dragStartY;

                int minAbove = getTrackMinimumHeight(abovePane);
                int newAboveHeight = Math.max(minAbove, originalAboveHeight + delta);
                setTrackHeight(abovePane, newAboveHeight);
            }
        };

        addMouseListener(mouseAdapter);
        addMouseMotionListener(mouseAdapter);
    }

    private static int getTrackHeight(TrackPanelScrollPane pane) {
        Track track = pane.getTrackPanel().getTrack();
        return track != null ? track.getHeight() : 0;
    }

    private static int getTrackMinimumHeight(TrackPanelScrollPane pane) {
        Track track = pane.getTrackPanel().getTrack();
        return track != null ? track.getMinimumHeight() : 10;
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

