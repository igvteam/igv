/*
* To change this template, choose Tools | Templates
* and open the template in the editor.
*/
package org.igv.ui.legend;

import org.igv.prefs.PreferencesManager;
import org.igv.track.TrackType;
import org.igv.ui.IGV;
import org.igv.ui.UIConstants;
import org.igv.ui.WaitCursorManager;
import org.igv.ui.util.IGVMouseInputAdapter;
import org.igv.ui.util.UIUtilities;

import javax.swing.*;
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

        mouseListener = new IGVMouseInputAdapter() {

            @Override
            public void mouseEntered(MouseEvent e) {
                LegendPanel.this.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
            }

            @Override
            public void mouseExited(MouseEvent e) {
            }

            @Override
            public void igvMouseClicked(MouseEvent e) {
                edit();
            }
        };
        addMouseListener(mouseListener);

        UIUtilities.invokeOnEventThread(() -> LegendPanel.this.setToolTipText(UIConstants.CLICK_ITEM_TO_EDIT_TOOLTIP));
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
        if (PreferencesManager.getPreferences().getAntiAliasing()) {
            ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        }
        paintLegend((Graphics2D)g);

    }

    abstract protected void paintLegend(Graphics2D g);

    abstract public void edit();

    abstract protected void reloadPreferences();

    abstract protected void persistResetPreferences();
}
