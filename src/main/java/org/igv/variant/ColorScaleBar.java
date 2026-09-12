package org.igv.variant;

import org.igv.renderer.ContinuousColorScale;
import org.igv.ui.FontManager;

import javax.swing.JPanel;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;

/**
 * A color scale drawn as a gradient with its values labelled -- the minimum and maximum, and for a double
 * gradient the edges of the neutral band.  Without the values a gradient says which colors are used but not
 * what any of them mean.
 */
public class ColorScaleBar extends JPanel {

    private static final int BAR_HEIGHT = 18;
    private static final int LABEL_HEIGHT = 14;

    private ContinuousColorScale scale;

    public ColorScaleBar(ContinuousColorScale scale, int width) {
        this.scale = scale;
        setPreferredSize(new Dimension(width, BAR_HEIGHT + LABEL_HEIGHT));
    }

    public void setScale(ContinuousColorScale scale) {
        this.scale = scale;
        repaint();
    }

    @Override
    protected void paintComponent(Graphics g) {

        super.paintComponent(g);
        if (scale == null) {
            return;
        }

        Graphics2D g2d = (Graphics2D) g;
        g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

        final double min = scale.getMinimum();
        final double max = scale.getMaximum();
        final int width = getWidth();

        for (int x = 0; x < width; x++) {
            double value = min + (max - min) * x / Math.max(1, width - 1);
            g2d.setColor(scale.getColor((float) value));
            g2d.drawLine(x, 0, x, BAR_HEIGHT);
        }

        g2d.setFont(FontManager.getFont(10));
        g2d.setColor(getForeground());
        FontMetrics metrics = g2d.getFontMetrics();
        int baseline = BAR_HEIGHT + metrics.getAscent() + 1;

        String minLabel = format(min, max - min);
        String maxLabel = format(max, max - min);
        g2d.drawString(minLabel, 0, baseline);
        g2d.drawString(maxLabel, width - metrics.stringWidth(maxLabel), baseline);

        if (scale.isUseDoubleGradient()) {
            int minEnd = metrics.stringWidth(minLabel);
            int maxStart = width - metrics.stringWidth(maxLabel);
            // The band edges coincide in the common case, so draw one label rather than two identical ones
            if (scale.getNegStart() == scale.getPosStart()) {
                drawInteriorLabel(g2d, metrics, scale.getNegStart(), min, max, baseline, minEnd, maxStart);
            } else {
                drawInteriorLabel(g2d, metrics, scale.getNegStart(), min, max, baseline, minEnd, maxStart);
                drawInteriorLabel(g2d, metrics, scale.getPosStart(), min, max, baseline, minEnd, maxStart);
            }
        }
    }

    /**
     * Draw a label at its position on the scale, unless it would collide with the end labels.
     */
    private void drawInteriorLabel(Graphics2D g2d, FontMetrics metrics, double value, double min, double max,
                                   int baseline, int minEnd, int maxStart) {

        String label = format(value, max - min);
        int centre = (int) ((value - min) / (max - min) * getWidth());
        int x = centre - metrics.stringWidth(label) / 2;

        if (x > minEnd + 4 && x + metrics.stringWidth(label) < maxStart - 4) {
            g2d.drawString(label, x, baseline);
            g2d.setColor(getForeground());
            g2d.drawLine(centre, BAR_HEIGHT, centre, BAR_HEIGHT + 2);
        }
    }

    /**
     * Show enough decimal places for the labels to differ.  Allele frequencies and phred scores differ by
     * orders of magnitude, so a fixed format would print "0.0" at both ends of one scale.
     */
    static String format(double value, double range) {
        if (range >= 100) {
            return String.format("%.0f", value);
        } else if (range >= 10) {
            return String.format("%.1f", value);
        } else if (range >= 1) {
            return String.format("%.2f", value);
        } else {
            return String.format("%.3f", value);
        }
    }
}
