package org.igv.sashimi;

import org.igv.ui.panel.Paintable;

import javax.swing.*;
import java.awt.*;

public class SashimiContentPane extends JSplitPane implements Paintable {

    JPanel sashimiPanel;
    JScrollPane scrollableGenePane;

    public SashimiContentPane(JPanel sashimiPanel, JScrollPane scrollableGenePane) {
        super(JSplitPane.VERTICAL_SPLIT, sashimiPanel, scrollableGenePane);
        this.sashimiPanel = sashimiPanel;
        this.scrollableGenePane = scrollableGenePane;
    }

    @Override
    public void paintOffscreen(Graphics2D g, Rectangle rect, boolean batch) {
        this.sashimiPanel.paint(g);
        g.translate(0, this.sashimiPanel.getHeight() + 10);
        this.scrollableGenePane.getViewport().paint(g);
    }

    @Override
    public int getSnapshotHeight(boolean batch) {
        return this.getHeight();
    }
}
