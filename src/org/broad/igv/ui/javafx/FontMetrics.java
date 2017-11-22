/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2017 Broad Institute
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

package org.broad.igv.ui.javafx;

import javafx.scene.text.Font;
import javafx.scene.text.Text;

/**
 * A substitute for the old Swing FontMetrics class since JavaFX doesn't have a (public) version.
 * This only implements the method(s) we need, so it's not general purpose.
 *
 * @author eby
 */
public class FontMetrics {

    // From https://stackoverflow.com/questions/13015698/how-to-calculate-the-pixel-width-of-a-string-in-javafx

    /**
     * Get the width of a text String with the given Font applied.
     * This is taken from https://stackoverflow.com/questions/13015698/how-to-calculate-the-pixel-width-of-a-string-in-javafx
     * Note: for now, this code takes the simple form of working with the font ONLY.  It may need to be more complex for
     * dealing with CSS.
     *
     * @param textString
     * @param font
     * @return
     */
    public static double getTextWidthInFont(String textString, Font font) {
        // We might consider caching these Text objects if it's too expensive to create them. 
        Text text = new Text(textString);
        text.setFont(font);
        return text.getBoundsInLocal().getWidth();
    }

    /**
     * Get the width of a text String with the Font of the given Text instance.  This is an optimization for where it's
     * necessary to calculate widths for many Strings.
     *
     * @param textString
     * @param text
     * @return
     */
    public static double getTextWidthInFont(String textString, Text text) {
        text.setText(textString);
        return text.getBoundsInLocal().getWidth();
    }

}
