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

import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.renderer.ColorScale;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.util.IGVMouseInputAdapter;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import javax.swing.event.MouseInputListener;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.function.Consumer;

/**
 * @author jrobinso
 */
abstract public class LegendPanel extends JPanel {

    protected final Set<Consumer<ColorScale>> changeListeners = new LinkedHashSet<>();

    /**
     * Constructs ...
     */
    public LegendPanel() {

        addMouseListener(new IGVMouseInputAdapter() {
            @Override
            public void mouseEntered(MouseEvent e) {
                LegendPanel.this.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
            }

            @Override
            public void igvMouseClicked(MouseEvent e) {
                edit();
            }
        });

        UIUtilities.invokeOnEventThread(() -> LegendPanel.this.setToolTipText(UIConstants.CLICK_ITEM_TO_EDIT_TOOLTIP));
    }

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
     * @param g
     */
    @Override
    public void paintComponent(Graphics g) {
        super.paintComponent(g);
      //  if (PreferencesManager.getPreferences().getAntiAliasing()) {
            ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        //}
        paintLegend((Graphics2D)g);

    }

    public void addChangeListener(Consumer<ColorScale> listener) {
        changeListeners.add(listener);
    }

    abstract protected void paintLegend(Graphics2D g);

    /**
     * Open the user preferences dialog
     */
    abstract public void edit();

    /**
     * reload preferences from the saved IGV preferences
     */
    abstract protected void reloadPreferences();

    /**
     * reset the IGV preferences to their default
     */
    abstract protected void resetPreferencesToDefault();

    /**
     * write back the current preferences in this dialog to the IGV preferences
     */
    abstract protected void persistCurrentPreferences();

    abstract public ColorScale getColorScale();
}
