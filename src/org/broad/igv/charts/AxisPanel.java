package org.broad.igv.charts;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.AffineTransform;
import java.text.DecimalFormat;

/**
 * @author Jim Robinson
 * @date 10/30/11
 */
class AxisPanel extends JComponent {

    //

    private static final int TICK_SIZE = 4;
    private static final int TICK_GAP = 2;
    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat();


    Axis axis;


    public void setAxisModel(Axis axis) {
        this.axis = axis;
    }


    @Override
    public void setBounds(int x, int y, int w, int h) {
        super.setBounds(x, y, w, h);
        updateAxisDimension(w, h);
    }

    private void updateAxisDimension(int w, int h) {

        if (axis == null) return;

        if (axis.getOrientation() == Axis.Orientation.HORIZONTAL) {
            axis.setPanelSize(w);
        } else {
            axis.setPanelSize(h);
        }
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);

        ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        if (axis == null) return;

        final Graphics2D g2D = (Graphics2D) g;
        final FontMetrics fontMetrics = g.getFontMetrics();
        final int minusStrWidth = fontMetrics.stringWidth("-") / 2;
        final int strHeight = fontMetrics.getHeight();
        final int top = 0;
        final int bottom = getHeight();
        final int width = getWidth();
        final double[] ticks = axis.ticks;
        final double tickStart = ticks[0];
        final double tickIncrement = ticks[1];

        double labelStep = 1;


        if (axis.getOrientation() == Axis.Orientation.HORIZONTAL) {
            int px;
            final int tickLabelY = TICK_SIZE + 2 * TICK_GAP + strHeight;
            double xtick = tickStart;
            do {
                px = axis.getPixelForValue(xtick);
                if (px > 0 && px < width) {
                    g.drawLine(px, TICK_GAP, px, TICK_GAP + TICK_SIZE);

                    //if (xtick % labelStep == 0) {
                    String label = DECIMAL_FORMAT.format(xtick);
                    if (label.equals("-0")) label = "0";
                    int strWidth = fontMetrics.stringWidth(label);
                    if (px > strWidth && (px + strWidth) < width) {
                        int strPosition = px - strWidth / 2 + 1;
                        if (xtick < 0) {
                            strPosition -= minusStrWidth;
                        }
                        g.drawString(label, strPosition, tickLabelY);
                    }
                    //}

                }
                xtick += tickIncrement;
            } while (px < width);

            String label = axis.getLabel();
            if (label != null) {
                int strWidth = fontMetrics.stringWidth(label);
                int strX = (getWidth() - strWidth) / 2;
                g.drawString(label, strX, bottom - 20);

            }
        } else {
            double ytick = tickStart;
            int py;
            do {
                py = axis.getPixelForValue(ytick);
                if (py - strHeight > top && py < bottom) {
                    g.drawLine(width - (TICK_GAP + TICK_SIZE), py, width - TICK_GAP, py);

                    //if (ytick % labelStep == 0) {
                    String label = DECIMAL_FORMAT.format(ytick);
                    if (label.equals("-0")) label = "0";
                    int strWidth = fontMetrics.stringWidth(label);
                    int strPosition = width - (TICK_GAP + TICK_SIZE) - 2 * TICK_GAP - strWidth;
                    int tickLabelY = py + strHeight / 2 - 3;
                    g.drawString(label, strPosition, tickLabelY);
                    //}


                }
                ytick += tickIncrement;
            } while (py > top);


            String label = axis.getLabel();
            if (label != null) {
                int strWidth = fontMetrics.stringWidth(label);
                int strX = (bottom + strWidth) / 2;

                AffineTransform t = ((Graphics2D) g).getTransform();
                AffineTransform rotateTransform = new AffineTransform();
                rotateTransform.quadrantRotate(-1);
                ((Graphics2D) g).transform(rotateTransform);
                g.drawString(label, -strX, 20);
                ((Graphics2D) g).setTransform(t);

            }

        }
    }
}
