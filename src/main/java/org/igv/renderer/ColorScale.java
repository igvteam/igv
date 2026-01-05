package org.igv.renderer;

import java.awt.*;

/**
 * @author jrobinso
 */
public interface ColorScale {

    public Color getColor(String symbol);

    public Color getColor(float value);

    public Color getNoDataColor();

    public String asString();

    public boolean isDefault();

}
