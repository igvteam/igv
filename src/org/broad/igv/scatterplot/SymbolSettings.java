package org.broad.igv.scatterplot;

import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Apr 6, 2010
 * Time: 5:05:05 PM
 * To change this template use File | Settings | File Templates.
 */
/*
*   Container class for plotting symbols assigned to ScatterPlotFrame series
*
*   Note: This class contains all information needed to plot a symbol
*       for the particular series it is assigned to.
* */
public class SymbolSettings {

    private int series;             // plot series the symbol settings are assigned to
    private String category;        // category
    private Shape plotShape;        // plot shape
    private Color plotColor;        // shape fill color
    private boolean isOutlined;     // boolean indicates if shape is outlined in black
    private boolean isFilled;       // boolean indicates shape is filled with plot color
    private boolean isVisible;      // plot shape as is visible?

    public SymbolSettings(int series, String category, Shape shape, Color plotColor,
                          boolean isOutlined, boolean isFilled, boolean isVisible){

            this. series = series;
            this.category = category;
            this.plotShape = shape;
            this.plotColor = plotColor;
            this.isOutlined = isOutlined;
            this.isFilled = isFilled;
            this.isVisible = isVisible;
    }

    public void setSeries(int series) {
        this.series = series;
    }

    public int getSeries() {
        return series;
    }

    public void setCategory(String category) {
       this.category = category;
    }

    public String getCategory() {
       return category;
    }

    public void setPlotShape(Shape plotShape) {
       this.plotShape = plotShape;
    }

    public Shape getPlotShape() {
       return plotShape;
    }

    public void setPlotColor(Color plotColor) {
       this.plotColor = plotColor;
    }

    public Color getPlotColor() {
       return plotColor;
    }

    public boolean isOutlined() {
           return isOutlined;
       }

    public void setOutlined(boolean outlined) {
        isOutlined = outlined;
    }

    public void setFilled(boolean isFilled) {
        this.isFilled = isFilled;
    }

    public boolean isFilled() {
        return isFilled;
    }

    public void setVisible(boolean isVisible) {
        this.isVisible = isVisible;
    }

    public boolean isVisible() {
       return isVisible;
    }

}
