package org.broad.igv.charts;

/**
 * @author Jim Robinson
 * @date 10/31/11
 */
public class XYDataPoint {

    double x;
    double y;
    String description;

    public XYDataPoint(double x, double y, String description) {
        this.x = x;
        this.y = y;
        this.description = description;
    }


    public double getX() {
        return x;
    }

    public double getY() {
        return y;
    }

    public String getDescription() {
        return description;
    }

    /**
     * Test if the point represented by (px, py) is with tolerance of this point.
     *
     *
     * @param px
     * @param py
     * @param toleranceX
     * @param toleranceY
     * @return
     */
    public boolean contains(double px, double py, double toleranceX, double toleranceY) {

        return px > x - toleranceX && px < x + toleranceX && py > y - toleranceY && py < y + toleranceY;

    }
}
