/*
 * AlphaColorGradient.java
 *
 * Created on November 17, 2007, 3:16 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package org.broad.igv.util;

import java.awt.*;

/**
 * @author jrobinso
 */
public class AlphaColorGradient {
    Color maxColor = Color.red;
    Color minColor = Color.blue;
    double minValue = -1.5;
    double midValue = 0;
    double maxValue = 1.5;

    Color getColor(double value) {
        if (value > midValue) {
            int alpha = (int) (255 * (value - midValue) / (maxValue - midValue));
            alpha = Math.max(0, Math.min(255, alpha));
            return new Color(maxColor.getRed(), maxColor.getGreen(), maxColor.getBlue(), alpha);
        } else {
            int alpha = (int) (255 * (midValue - value) / (midValue - minValue));
            alpha = Math.max(0, Math.min(255, alpha));
            return new Color(minColor.getRed(), minColor.getGreen(), minColor.getBlue(), alpha);
        }
    }
}

