/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */


/*
* To change this template, choose Tools | Templates
* and open the template in the editor.
*/
package org.broad.igv.renderer;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.color.ColorUtilities;

import java.awt.*;
import java.awt.image.BufferedImage;

/**
 * @author jrobinso
 */
public class ContinuousColorScale extends AbstractColorScale {


    /**
     * String to use to identify this class when serializing.  We could use
     * the class name, but that would invalidate serialized instances if the
     * name was ever changed.
     */
    public static final String serializedClassName = "ContinuousColorScale";
    private final boolean useDoubleGradient;
    private double negEnd;
    private double posEnd;
    private double negStart;
    private double posStart;
    private Color minColor;
    private Color midColor = Color.white;
    private Color maxColor;
    private Color[] colors;
    private double delta;

    /**
     * Constructs ...
     *
     * @param string
     */

    public ContinuousColorScale(String string) {

        String[] tokens = string.split(";");
        if (tokens.length == 5) {
            this.negEnd = Double.parseDouble(tokens[1]);
            this.posEnd = Double.parseDouble(tokens[2]);
            this.minColor = ColorUtilities.stringToColor(tokens[3]);
            this.maxColor = ColorUtilities.stringToColor(tokens[4]);
            this.useDoubleGradient = false;
        } else if (tokens.length == 8) {
            this.negStart = Double.parseDouble(tokens[1]);
            this.negEnd = Double.parseDouble(tokens[2]);
            this.posStart = Double.parseDouble(tokens[3]);
            this.posEnd = Double.parseDouble(tokens[4]);
            this.minColor = ColorUtilities.stringToColor(tokens[5]);
            this.midColor = ColorUtilities.stringToColor(tokens[6]);
            this.maxColor = ColorUtilities.stringToColor(tokens[7]);
            this.useDoubleGradient = true;
        } else {
            throw new RuntimeException("Illegal ColorScale: " + string);
        }
        initColors();


    }


    /**
     * Constructs ...
     *
     * @param min
     * @param max
     * @param minColor
     * @param maxColor
     */
    public ContinuousColorScale(double min, double max, Color minColor, Color maxColor) {
        this.negEnd = min;
        this.posEnd = max;
        this.negStart = Math.max(0, min);
        this.posStart = this.negStart;
        this.minColor = minColor;
        this.maxColor = maxColor;
        this.useDoubleGradient = false;
        initColors();
    }

    /**
     * Constructs ...
     *
     * @param min
     * @param mid
     * @param max
     * @param minColor
     * @param midColor
     * @param maxColor
     */
    public ContinuousColorScale(double min, double mid, double max, Color minColor, Color midColor,
                                Color maxColor) {
        this.negEnd = min;
        this.posEnd = max;
        this.negStart = mid;
        this.posStart = mid;
        this.minColor = minColor;
        this.midColor = midColor;
        this.maxColor = maxColor;
        this.useDoubleGradient = true;
    }

    /**
     * Constructs ...
     *
     * @param negStart
     * @param negEnd
     * @param posStart
     * @param posEnd
     * @param minColor
     * @param midColor
     * @param maxColor
     */
    public ContinuousColorScale(double negStart, double negEnd, double posStart, double posEnd,
                                Color minColor, Color midColor, Color maxColor) {
        this.negEnd = negEnd;
        this.posEnd = posEnd;
        this.negStart = negStart;
        this.posStart = posStart;
        this.minColor = minColor;
        this.midColor = midColor;
        this.maxColor = maxColor;
        this.useDoubleGradient = true;

    }

    public ContinuousColorScale(ContinuousColorScale otherScale) {
        this.negEnd = otherScale.negEnd;
        this.posEnd = otherScale.posEnd;
        this.negStart = otherScale.negStart;
        this.posStart = otherScale.posStart;
        this.minColor = otherScale.minColor;
        this.midColor = otherScale.midColor;
        this.maxColor = otherScale.maxColor;
        this.useDoubleGradient = true;
    }

    public void setNegEnd(double negEnd) {
        this.negEnd = negEnd;

    }

    public void setPosEnd(double posEnd) {
        this.posEnd = posEnd;
        colors = null;
    }

    public void setNegStart(double negStart) {
        this.negStart = negStart;
        colors = null;
    }

    public void setPosStart(double posStart) {
        this.posStart = posStart;
        colors = null;
    }

    public void setMinColor(Color minColor) {
        this.minColor = minColor;
        colors = null;
    }

    public void setMidColor(Color midColor) {
        this.midColor = midColor;
        colors = null;
    }

    public void setMaxColor(Color maxColor) {
        this.maxColor = maxColor;
        colors = null;
    }

    /**
     * Create a string form of this object.
     *
     * @return
     */
    public String asString() {
        StringBuilder buf = new StringBuilder();
        buf.append(serializedClassName + ";");
        if (useDoubleGradient) {
            buf.append(String.valueOf(negStart) + ";");
            buf.append(String.valueOf(negEnd) + ";");
            buf.append(String.valueOf(posStart) + ";");
            buf.append(String.valueOf(posEnd) + ";");
            buf.append(ColorUtilities.colorToString(minColor) + ";");
            buf.append(ColorUtilities.colorToString(midColor) + ";");
            buf.append(ColorUtilities.colorToString(maxColor));
        } else {
            buf.append(String.valueOf(negEnd) + ";");
            buf.append(String.valueOf(posEnd) + ";");
            buf.append(ColorUtilities.colorToString(minColor) + ";");
            buf.append(ColorUtilities.colorToString(maxColor));
        }
        return buf.toString();
    }



    private void initColors() {
        colors = new Color[251];
        delta = (posEnd - negEnd) / colors.length;
        if (isUseDoubleGradient()) {
            ColorGradient csPos = new ColorGradient(posStart, posEnd, midColor, maxColor);
            ColorGradient csNeg = new ColorGradient(negEnd, negStart, minColor, midColor);

            for (int i = 0; i < colors.length; i++) {
                double x = getMinimum() + i * delta;
                if ((x > negStart) && (x < posStart)) {
                    colors[i] = midColor;
                } else if (x <= negStart) {
                    colors[i] = csNeg.getColor(x);
                } else {
                    colors[i] = csPos.getColor(x);
                }
            }

        } else {

            ColorGradient cs = new ColorGradient(negEnd, posEnd, minColor, maxColor);
            for (int i = 0; i < colors.length; i++) {
                double x = getMinimum() + i * delta;
                colors[i] = cs.getColor(x);
            }
        }
    }

    /**
     * Method description
     *
     * @param val
     * @return
     */
    @Override
    public Color getColor(float val) {

        if (colors == null) {
            initColors();
        }

        if(Float.isNaN(val)) {
            return PreferencesManager.getPreferences().getAsColor(Constants.NO_DATA_COLOR);
        }

        // See if we are in the midrange.

        if (val >= negStart && val <= posStart) {
            return midColor;
        } else {
            //double f = (val - getMinimum()) / (getMaximum() - getMinimum());
            int index = (int) Math.round((val - negEnd) / delta);
            index = Math.max(0, Math.min(index, colors.length - 1));
            return colors[index];
        }
    }

    /**
     * Method description
     *
     * @return
     */
    public double getNegStart() {
        return negStart;
    }

    /**
     * Method description
     *
     * @return
     */
    public double getPosStart() {
        return posStart;
    }

    /**
     * Method description
     *
     * @return
     */
    public double getMinimum() {
        return negEnd;
    }

    /**
     * Method description
     *
     * @return
     */
    public double getMaximum() {
        return posEnd;
    }

    public double getBaseline() {
        return (negStart + posStart) / 2;
    }

    /**
     * Method description
     *
     * @return
     */
    public Color getMinColor() {
        return minColor;
    }

    /**
     * Method description
     *
     * @return
     */
    public Color getMidColor() {
        return midColor;
    }

    /**
     * Method description
     *
     * @return
     */
    public Color getMaxColor() {
        return maxColor;
    }

    /**
     * Method description
     *
     * @return
     */
    public boolean isUseDoubleGradient() {
        return useDoubleGradient;
    }

    static class ColorGradient {

        private final Color minColor;
        private final Color maxColor;

        private final double min;
        private final double max;
        private final BufferedImage posImage;


        /**
         * Constructs ...
         *
         * @param min
         * @param max
         * @param negColor
         * @param posColor
         */
        public ColorGradient(double min, double max, Color negColor, Color posColor) {
            this.min = min;
            this.max = max;
            this.minColor = negColor;
            this.maxColor = posColor;
            this.posImage = createGradientImage(negColor, posColor);
        }

        /**
         * Creates a gradient image given specified <CODE>Color</CODE>(s)
         *
         * @param color1 <CODE>Color</CODE> to display at left side of gradient
         * @param color2 <CODE>Color</CODE> to display at right side of gradient
         * @return returns a gradient image
         */
        private static BufferedImage createGradientImage(Color color1, Color color2) {

            BufferedImage image =
                    GraphicsEnvironment.getLocalGraphicsEnvironment()
                            .getDefaultScreenDevice().getDefaultConfiguration()
                            .createCompatibleImage(256, 1);
            Graphics2D graphics = image.createGraphics();
            GradientPaint gp = new GradientPaint(0, 0, color1, 255, 0, color2);

            graphics.setPaint(gp);
            graphics.drawRect(0, 0, 255, 1);
            graphics.dispose();

            return image;
        }

        public Color getColor(double value) {
            if(value <= this.min) {
                return minColor;
            } else if(value >= this.max) {
                return maxColor;
            }

            double span = this.max - this.min;
            int colorIndex = 0;

            if (value <= min) {
                colorIndex = 0;
            } else if (value >= max) {
                colorIndex = 255;
            } else {
                colorIndex = (int) (((value - this.min) / span) * 255);
            }

            int rgb = posImage.getRGB(colorIndex, 0);

            return new Color(rgb);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ContinuousColorScale that = (ContinuousColorScale) o;

        if (Double.compare(that.negEnd, negEnd) != 0) return false;
        if (Double.compare(that.negStart, negStart) != 0) return false;
        if (Double.compare(that.posEnd, posEnd) != 0) return false;
        if (Double.compare(that.posStart, posStart) != 0) return false;
        if (useDoubleGradient != that.useDoubleGradient) return false;
        if (!maxColor.equals(that.maxColor)) return false;
        if (!midColor.equals(that.midColor)) return false;
        if (!minColor.equals(that.minColor)) return false;
        if (!noDataColor.equals(that.noDataColor)) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = (useDoubleGradient ? 1 : 0);
        temp = negEnd != +0.0d ? Double.doubleToLongBits(negEnd) : 0L;
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = posEnd != +0.0d ? Double.doubleToLongBits(posEnd) : 0L;
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = negStart != +0.0d ? Double.doubleToLongBits(negStart) : 0L;
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = posStart != +0.0d ? Double.doubleToLongBits(posStart) : 0L;
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        result = 31 * result + minColor.hashCode();
        result = 31 * result + midColor.hashCode();
        result = 31 * result + maxColor.hashCode();
        result = 31 * result + noDataColor.hashCode();
        return result;
    }
}
