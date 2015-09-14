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
package org.broad.igv.ui.legend;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.PreferenceManager;
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
        if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.ENABLE_ANTIALISING)) {
            ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        }
        paintLegend(g);
        ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_DEFAULT);

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
