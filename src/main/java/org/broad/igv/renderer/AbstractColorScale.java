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

import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;

import java.awt.*;

/**
 * @author jrobinso
 */
public abstract class AbstractColorScale implements ColorScale {

    final static Color DEFAULT_COLOR = Color.BLACK;
    protected Color noDataColor = PreferencesManager.getPreferences().getAsColor(Constants.NO_DATA_COLOR);

    /* default behavior is parse the input as a float and look up the color by number*/
    public Color getColor(String symbol) {
        try{
            return getColor(Float.parseFloat(symbol));
        } catch (NumberFormatException |  NullPointerException e){
            return DEFAULT_COLOR;
        }
    }

    public Color getColor(float value) {
        return DEFAULT_COLOR;
    }

    /**
     * Method description
     *
     * @param color
     */
    public void setNoDataColor(Color color) {
        this.noDataColor = color;

    }

    public Color getNoDataColor() {
        return noDataColor;
    }
}
