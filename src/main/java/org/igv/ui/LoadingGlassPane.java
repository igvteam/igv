package org.igv.ui;

import javax.swing.*;
import java.awt.*;

/**
 * The main window's glass pane.  Shows a message in the center of the window while a long running task that
 * disables the menus, such as loading a genome, is in progress, so it is apparent why the menus are unavailable.
 * The message appears only after a short delay, so a task that finishes quickly does not flash it.  The glass pane
 * is made visible by the WaitCursorManager; with no message it paints nothing.
 */
public class LoadingGlassPane extends JComponent {

    private static final int SHOW_DELAY_MS = 500;

    private String message;
    private boolean showing;   // true once the delay has passed
    private final Timer showTimer;

    public LoadingGlassPane() {
        showTimer = new Timer(SHOW_DELAY_MS, e -> {
            showing = true;
            repaint();
        });
        showTimer.setRepeats(false);
    }

    /**
     * @param message the message to show, or null to remove it
     */
    public void setMessage(String message) {
        this.message = message;
        if (message == null) {
            showTimer.stop();
            showing = false;
        } else if (!showing && !showTimer.isRunning()) {
            showTimer.start();
        }
        repaint();
    }

    @Override
    protected void paintComponent(Graphics g) {
        if (message == null || !showing) {
            return;
        }
        Graphics2D g2 = (Graphics2D) g.create();
        try {
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            g2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
            g2.setFont(UIManager.getFont("Label.font").deriveFont(Font.PLAIN, 14f));
            FontMetrics fm = g2.getFontMetrics();

            // A box around the message, in the look and feel's own colors
            int w = Math.min(getWidth() - 20, fm.stringWidth(message) + 30);
            int h = fm.getHeight() + 20;
            int x = (getWidth() - w) / 2;
            int y = (getHeight() - h) / 2;
            g2.setColor(UIManager.getColor("Panel.background"));
            g2.fillRoundRect(x, y, w, h, 12, 12);
            Color textColor = UIManager.getColor("TitlePane.embeddedForeground");
            g2.setColor(textColor != null ? textColor : UIManager.getColor("Label.foreground"));
            g2.drawRoundRect(x, y, w, h, 12, 12);

            g2.clipRect(x, y, w, h);
            g2.drawString(message, x + (w - fm.stringWidth(message)) / 2, y + (h - fm.getHeight()) / 2 + fm.getAscent());
        } finally {
            g2.dispose();
        }
    }
}
