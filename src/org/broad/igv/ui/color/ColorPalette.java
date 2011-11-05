package org.broad.igv.ui.color;

import java.awt.*;

/**
 * @author Jim Robinson
 * @date 11/5/11
 */
public class ColorPalette  {

    private String name;
    private Color[]  colors;

    public ColorPalette(String name, Color[] colors) {
        this.colors = colors;
        this.name = name;
    }

    public String getName() {
        return name;
    }

    public Color[] getColors() {
        return colors;
    }
}
