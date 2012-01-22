package org.broad.igv.charts;

/**
 * @author Jim Robinson
 * @date 10/31/11
 */
public class XYDataPoint {

    double x;
    double y;
    int mutationCount;   // <= TODO make this some generic attribute, we've lost the generality of this class
    private double methylation = Double.NaN; // DITTO
    private double expression = Double.NaN;  // DITTO
    private double copyNumber = Double.NaN;  // DITTO
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

    public int getMutationCount() {
        return mutationCount;
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

    public void setMutationCount(int mutationCount) {
        this.mutationCount = mutationCount;
    }

    public void setCopyNumber(double copyNumber) {
        this.copyNumber = copyNumber;
    }


    public double getMethylation() {
        return methylation;
    }

    public void setMethylation(double methylation) {
        this.methylation = methylation;
    }

    public double getExpression() {
        return expression;
    }

    public void setExpression(double expression) {
        this.expression = expression;
    }

    public double getCopyNumber() {
        return copyNumber;
    }
}
