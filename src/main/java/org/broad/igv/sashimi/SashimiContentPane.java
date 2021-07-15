package org.broad.igv.sashimi;

import org.broad.igv.ui.panel.Paintable;

import javax.swing.*;
import java.awt.*;

public class SashimiContentPane extends JPanel implements Paintable {

    @Override
    public void paintOffscreen(Graphics2D g, Rectangle rect, boolean batch) {
        this.paint(g);
    }

    @Override
    public int getSnapshotHeight(boolean batch) {
        return this.getHeight();
    }
}
