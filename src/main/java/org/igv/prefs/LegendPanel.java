package org.igv.prefs;

import org.igv.Globals;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.renderer.ContinuousColorScale;
import org.igv.track.TrackType;
import org.igv.ui.FontManager;
import org.igv.ui.IGV;
import org.igv.ui.UIConstants;
import org.igv.ui.legend.HeatmapLegendEditor;
import org.igv.ui.util.IGVMouseInputAdapter;
import org.igv.ui.util.UIUtilities;

import javax.swing.*;
import javax.swing.event.MouseInputListener;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.text.DecimalFormat;

/**
 * @author jrobinso
 * A panel that displays a legend for a heatmap track in the preferences editor.  Clicking on the panel
 * brings up the heatmap editor dialog.
 */
public class LegendPanel extends JPanel {

    static Logger log = LogManager.getLogger(LegendPanel.class);

    enum Orientation {HORIZONTAL, VERTICAL}

    private Orientation orientation = Orientation.HORIZONTAL;
    private String key;
    protected ContinuousColorScale colorScale;
    private MouseInputListener mouseListener;

    public LegendPanel(String key) {

        this.key = key;

        this.colorScale = PreferencesManager.getPreferences().getColorScale(key);

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

    public void edit() {

        UIUtilities.invokeOnEventThread(() -> {

            IGV.getInstance().setStatusBarMessage("Setting view properties...");

            HeatmapLegendEditor dialog = new HeatmapLegendEditor(IGV.getInstance().getMainFrame(), true, colorScale);

            dialog.setTitle("HeatMap Preferences");
            dialog.setVisible(true);

            if (dialog.isCanceled()) {
                IGV.getInstance().resetStatusMessage();
                return;
            }

            colorScale = dialog.getColorScheme();
            PreferencesManager.getPreferences().setColorScale(key, colorScale);

            LegendPanel.this.repaint();
            IGV.getInstance().repaint();

        });
    }

    @Override
    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        if (PreferencesManager.getPreferences().getAntiAliasing()) {
            ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        }
        paintLegend((Graphics2D) g);
    }

    protected void paintLegend(Graphics2D g) {
        if (orientation == Orientation.HORIZONTAL) {
            paintHorizontal(g);
        } else {
            paintVertical(g);
        }
    }

    protected void paintHorizontal(Graphics2D g2D) {

        DecimalFormat formatter = new DecimalFormat("0.0");

        g2D.setFont(FontManager.getFont(10));

        int npts = 5;
        double max = colorScale.getMaximum();
        double min = colorScale.getMinimum();

        int w = getWidth() - 20;
        double dx = ((double) w) / npts;
        double dxj = dx / 10;
        double delta = (max - min) / npts;
        double deltaj = delta / 10;

        for (int i = 0; i < npts + 1; i++) {
            for (int j = i * 10; j < i * 10 + 10; j++) {
                double val = min + j * deltaj;

                Color c = colorScale.getColor((float) val);

                g2D.setColor(c);

                int x0 = (int) (j * dxj);
                int x1 = (int) ((j + 1) * dxj);

                g2D.fillRect(x0, 0, (x1 - x0), (int) (getHeight() / 2));
            }

            double labelVal = min + i * delta;
            int x0 = (int) (i * dx);

            g2D.setColor(Globals.isDarkMode() ? Color.WHITE : Color.BLACK);
            g2D.drawString(formatter.format(labelVal), x0, (int) getHeight() - 5);
        }

    }

    void paintVertical(Graphics2D g2D) {

        DecimalFormat formatter = new DecimalFormat("0.0");

        g2D.setFont(FontManager.getFont(10));

        int npts = 5;
        double max = colorScale.getMaximum();
        double min = colorScale.getMinimum();

        int h = getWidth() - 20;
        double dy = ((double) h) / npts;
        double dyj = dy / 10;
        double delta = (max - min) / npts;
        double deltaj = delta / 10;

        int x0 = 10;
        int dx = 10;
        int y0;
        int y1 = 0;

        for (int i = 0; i < npts + 1; i++) {
            for (int j = i * 10; j < i * 10 + 10; j++) {
                double val = min + j * deltaj;

                Color c = colorScale.getColor((float) val);

                g2D.setColor(c);

                y0 = (int) (j * dyj);
                y1 = (int) ((j + 1) * dyj);

                g2D.fillRect(x0, y0, dx, y1 - y0);
            }

            double labelVal = min + i * delta;

            g2D.setColor(Color.BLACK);
            g2D.drawString(formatter.format(labelVal), x0 + 15, y1 - 5);
        }
    }

}
