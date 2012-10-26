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


/*
* To change this template, choose Tools | Templates
* and open the template in the editor.
*/
package org.broad.igv.ui.legend;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.track.TrackType;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import javax.swing.event.MouseInputAdapter;
import javax.swing.event.MouseInputListener;
import java.awt.*;
import java.awt.event.MouseEvent;

/**
 * @author jrobinso
 */
abstract public class LegendPanel extends JPanel {

    protected TrackType type;
    private MouseInputListener mouseListener;
    private WaitCursorManager.CursorToken token;

    /**
     * Constructs ...
     */
    public LegendPanel() {

        mouseListener = new MouseInputAdapter() {

            @Override
            public void mouseEntered(MouseEvent e) {

                LegendPanel.this.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
            }

            @Override
            public void mouseExited(MouseEvent e) {
            }

            @Override
            public void mouseClicked(MouseEvent e) {
                edit();
            }
        };
        addMouseListener(mouseListener);

        UIUtilities.invokeOnEventThread(new Runnable() {

            public void run() {
                LegendPanel.this.setToolTipText(UIConstants.CLICK_ITEM_TO_EDIT_TOOLTIP);
            }
        });
    }

    abstract protected void resetPreferencesToDefault();

    protected void showResetDisplay() {
        try {
            reloadPreferences();

        }
        finally {
            IGV.getInstance().resetStatusMessage();
        }

    }

    /**
     * Method description
     *
     * @param type
     */
    public void setTrackType(TrackType type) {
        this.type = type;
    }

    /**
     * Method description
     *
     * @param g
     */
    @Override
    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        paintLegend(g);
    }

    abstract protected void paintLegend(Graphics g);

    /**
     * Open the user preferences dialog
     */
    abstract public void edit();

    /*
     * protected void setColorScheme(double minimum, double median, double maximum,
     * Color minColor, Color medianColor, Color maxColor) {
     * colorScheme =
     * new ContinuousColorScale(minimum, maximum,
     * minColor, medianColor, maxColor, median, median);
     * }
     */

    // abstract protected LinkedHashMap<String, PreferenceDescriptor> addPreferences();

    abstract protected void reloadPreferences();

    abstract protected void persistResetPreferences();
}
