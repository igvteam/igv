/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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

package org.broad.igv.ui.util;

import javax.swing.*;
import java.util.HashMap;
import java.util.Map;

public class IconFactory {

    public enum IconID {

        DRAG_AND_DROP,
        SLIDER,
        SLIDER_BOTTOM,
        SLIDER_TOP,
        SLIDER_BAR,
        ZOOM,
        ZOOM_PLUS,
        ZOOM_MINUS,
        ZOOM_IN,
        ZOOM_OUT,
        REGION_OF_INTEREST,
        NO_TOOLTIP,
        TOOLTIP,
        SEARCH,
        OPEN_HAND,
        REFRESH,
        HOME,
        FIST
    }

    private Map<IconID, ImageIcon> icons;
    private static IconFactory theInstance = null;

    public static IconFactory getInstance() {
        if (theInstance == null) {
            theInstance = new IconFactory();
        }
        return theInstance;

    }

    private IconFactory() {

        icons = new HashMap();
        icons.put(IconID.DRAG_AND_DROP,
                createImageIcon("/images/dragNdrop.png", "drag and drop"));
        icons.put(IconID.ZOOM,
                createImageIcon("/images/zoomin.gif", "zoom"));
        icons.put(IconID.ZOOM_IN,
                createImageIcon("/toolbarButtonGraphics/general/ZoomIn24.gif", "Zoom in"));
        icons.put(IconID.ZOOM_OUT,
                createImageIcon("/toolbarButtonGraphics/general/ZoomOut24.gif", "zoom out"));
        icons.put(IconID.SLIDER,
                createImageIcon("/images/slider.png", "zoom out"));
        icons.put(IconID.SLIDER_BOTTOM,
                createImageIcon("/images/dsliderbottom.png", "zoom out"));
        icons.put(IconID.SLIDER_TOP,
                createImageIcon("/images/dslidertop.png", "zoom out"));
        icons.put(IconID.SLIDER_BAR,
                createImageIcon("/images/dsliderbar.png", "zoom out"));
        icons.put(IconID.ZOOM_PLUS,
                createImageIcon("/images/zoom-plus.png", "zoom out"));
        icons.put(IconID.ZOOM_MINUS,
                createImageIcon("/images/zoom-minus.png", "zoom out"));
        icons.put(IconID.OPEN_HAND,
                createImageIcon("/images/cursor-openhand-dot.png", "drag"));
        icons.put(IconID.FIST,
                createImageIcon("/images/cursor-fist.png", "drag"));
        icons.put(IconID.REGION_OF_INTEREST,
                createImageIcon("/images/regionTool.png", "region of interest"));
        icons.put(IconID.REFRESH,
                createImageIcon("/toolbarButtonGraphics/general/Refresh24.gif", "refresh"));
        icons.put(IconID.HOME,
                createImageIcon("/toolbarButtonGraphics/navigation/Home24.gif", "home"));
        icons.put(IconID.NO_TOOLTIP,
                createImageIcon("/images/no-tooltip.png", "suppress tooltip"));
        icons.put(IconID.TOOLTIP,
                createImageIcon("/images/tooltip.png", "tooltip"));
    }

    public ImageIcon getIcon(IconID id) {
        return icons.get(id);

    }

    /**
     * Returns an ImageIcon at the specified URL, or null if the url was invalid.
     */
    private ImageIcon createImageIcon(String url, String description) {
        java.net.URL imgURL = getClass().getResource(url);
        if (imgURL != null) {
            return new ImageIcon(imgURL, description);
        } else {
            System.err.println("Couldn't find file: " + url);
            return null;
        }
    }
}
