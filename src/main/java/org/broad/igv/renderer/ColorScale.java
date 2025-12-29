/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.renderer;

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
