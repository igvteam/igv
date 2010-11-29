/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.ui;

import java.awt.*;

/**
 * @author eflakes
 */
public class FontManager {

    // Constants
    final public static String SCALABLE_FONT_NAME = "Lucida Sans";

    static public void applyScalableTextFont(Graphics2D graphics2) {

        Font font = getScalableFont(graphics2.getFont().getSize());
        //applyScalableTextHints(graphics2);
        graphics2.setFont(font);
    }

    static public void applyScalableTextFont(Graphics2D graphics2, int size) {

        Font font = getScalableFont(size);
        //applyScalableTextHints(graphics2);
        graphics2.setFont(font);
    }

    static public Font getScalableFont(int size) {
        return new Font(SCALABLE_FONT_NAME, Font.PLAIN, size);
    }

    static public Font getScalableFont(int attribute, int size) {
        return new Font(SCALABLE_FONT_NAME, attribute, size);
    }

    static public void applyScalableTextHints(Graphics2D graphics2) {

        graphics2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                RenderingHints.VALUE_ANTIALIAS_ON);
    }
}
