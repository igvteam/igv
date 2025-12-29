package org.broad.igv.ui.dnd;

import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;

public class GhostGlassPane extends JPanel {
    private AlphaComposite composite;
    private BufferedImage dragged = null;
    private Point location = new Point(0, 0);

    public GhostGlassPane() {
        setOpaque(false);
        composite = AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.5f);
    }

    public void setImage(BufferedImage dragged) {
        this.dragged = dragged;
    }

    public void setPoint(Point location) {
        this.location = location;
    }

    public void paintComponent(Graphics g) {
        if (dragged == null)
            return;

        Graphics2D g2 = (Graphics2D) g;
        g2.setComposite(composite);
        g2.drawImage(dragged,
                (int) (location.getX() - (dragged.getWidth(this) / 2)),
                (int) (location.getY() - (dragged.getHeight(this) / 2)),
                null);
    }
}
