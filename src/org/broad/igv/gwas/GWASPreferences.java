package org.broad.igv.gwas;

import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: jussi
 * Date: Nov 30, 2009
 * Time: 4:49:46 PM
 * To change this template use File | Settings | File Templates.
 */
public class GWASPreferences {


    private int minPointSize = 2;
    private int maxPointSize = 7;
    private boolean chrColors = true;
    private boolean alternatingColors = false;
    private boolean log10 = true;
    private boolean ln = false;
    private Color primaryColor = new Color(70, 105, 131);
    private Color secondaryColor = new Color(199, 81, 39);


    public boolean isAlternatingColors() {
        return alternatingColors;
    }

    public void setAlternatingColors(boolean alternatingColors) {
        this.alternatingColors = alternatingColors;
    }


    public boolean isChrColors() {
        return chrColors;
    }

    public void setChrColors(boolean chrColors) {
        this.chrColors = chrColors;
    }

    public Color getPrimaryColor() {
        return primaryColor;
    }

    public void setPrimaryColor(Color primaryColor) {
        this.primaryColor = primaryColor;
    }

    public Color getSecondaryColor() {
        return secondaryColor;
    }

    public void setSecondaryColor(Color secondaryColor) {
        this.secondaryColor = secondaryColor;
    }


    public int getMinPointSize() {
        return minPointSize;
    }

    public void setMinPointSize(int minPointSize) {
        this.minPointSize = minPointSize;
    }

    public int getMaxPointSize() {
        return maxPointSize;
    }

    public void setMaxPointSize(int maxPointSize) {
        this.maxPointSize = maxPointSize;
    }
}
