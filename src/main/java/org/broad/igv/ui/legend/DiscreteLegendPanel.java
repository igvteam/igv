package org.broad.igv.ui.legend;

import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.ColorPalette;
import org.broad.igv.ui.color.ColorTable;
import org.broad.igv.ui.color.PaletteColorTable;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import java.awt.*;
import java.util.Map;

public class DiscreteLegendPanel extends LegendPanel {
    public static final int WIDTH = 10;
    public static final int SWATCH_SEP = 20;
    public static final int LABEL_SEP = 40;
    protected PaletteColorTable colorTable;

    public DiscreteLegendPanel(PaletteColorTable colorTable){
        this.colorTable = colorTable;
    }

    protected void reloadPreferences() {
        repaint();
    }

    protected void persistCurrentPreferences() {
    }

    @Override
    public ColorTable getColorScale() {
        return colorTable;
    }

    @Override
    protected void resetPreferencesToDefault() {
    }

    /**
     * Open the user preferences dialog
     */
    public void edit() {

        UIUtilities.invokeOnEventThread(() -> {

            IGV.getInstance().setStatusBarMessage("Setting view properties...");

            DiscreteColorScaleEditor dialog = new DiscreteColorScaleEditor(IGV.getInstance().getMainFrame(), true, colorTable);

            dialog.setTitle("Preferences");
            dialog.pack();
            dialog.setVisible(true);


            if (dialog.isCancelled()) {
                IGV.getInstance().resetStatusMessage();
                return;
            }

            colorTable = dialog.getColorScale();
            changeListeners.forEach(c -> c.accept(colorTable));
            //PreferencesManager.getPreferences().setColorScale(type, colorScale);
            IGV.getInstance().repaint();
            try {
                reloadPreferences();
            } finally {
                UIUtilities.invokeOnEventThread(() -> SwingUtilities.getWindowAncestor(DiscreteLegendPanel.this).toFront());
                IGV.getInstance().resetStatusMessage();
            }
        });
    }

    @Override
    public void paintLegend(Graphics2D g2D) {

        if (colorTable == null) {
            return;
        }

        g2D.setFont(FontManager.getFont(10));


        FontMetrics fm = g2D.getFontMetrics();
        int dh = fm.getHeight() / 2 + 3;

        int x = 0;
        int lineHeight = 12;
        int y = getHeight() < 2 * lineHeight ? 0 : lineHeight;
        int colCount = 0;

        int pairs = colorTable.entrySet().size();

        ColorPalette palette = colorTable.getPalette();
        int paletteLength = palette != null ? palette.colors().length : 0;
        int extras = Math.max(paletteLength - pairs, 0);
        for (Map.Entry<String, Color> entry : colorTable.entrySet()) {
            String label = entry.getKey().replace("_", " "); //is this necessary anymore?
            int labelWidth = (int) fm.getStringBounds(label, g2D).getWidth();
            boolean insufficientSpace = x + labelWidth + WIDTH + SWATCH_SEP > getWidth();
            if(insufficientSpace){
                boolean moreLinesAvailable = y + lineHeight + 5 + WIDTH < getHeight();
                if(moreLinesAvailable){
                    y += lineHeight + 5; // next line
                    x = 0;
                } else {
                    g2D.drawString("...", x, y + dh);
                    break;
                }
            } else {
                drawSwatch(g2D, x, y, entry.getValue());
                g2D.drawString(label, x + SWATCH_SEP, y + dh);
                x += labelWidth + LABEL_SEP;
            }
        }

        if(extras > 0) {
            for (Color color : palette.colors()) {
                boolean insufficientSpace = WIDTH + SWATCH_SEP > getWidth();
                if (insufficientSpace) {
                    boolean moreLinesAvailable = y + lineHeight + 5 + WIDTH < getHeight();
                    if (moreLinesAvailable) {
                        y += lineHeight + 5; // next line
                        x = 0;
                    } else {
                        g2D.drawString("...", x, y + dh);
                        break;
                    }
                } else {
                    drawSwatch(g2D, x, y, color);
                    x += SWATCH_SEP;
                }
            }
        }

    }

    private static void drawSwatch(Graphics2D g2D, int x, int y, Color color) {
        g2D.setColor(color);
        g2D.fillRect(x, y, WIDTH, WIDTH);
        g2D.setColor(Color.BLACK);
        g2D.drawRect(x, y, WIDTH, WIDTH);
    }
}
