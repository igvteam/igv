/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.igv.ui.legend;

import org.igv.track.TrackType;
import org.igv.ui.FontManager;

import java.awt.*;

/**
 * @author jrobinso
 */
public class LohLegendPanel extends HeatmapLegendPanel {


    public LohLegendPanel() {
        super(TrackType.LOH);
    }


    @Override
    public void paintLegend(Graphics2D g2D) {

        g2D.setFont(FontManager.getFont(10));
        FontMetrics fm = g2D.getFontMetrics();
        int dh = fm.getHeight() / 2 + 3;


        int x = 0;
        int y = getHeight() / 2;

        String label = "Loss";
        int labelWidth = (int) fm.getStringBounds(label, g2D).getWidth();
        g2D.setColor(colorScale.getMaxColor());
        g2D.fillRect(x, y, 10, 10);
        g2D.setColor(Color.BLACK);
        g2D.drawRect(x, y, 10, 10);
        g2D.drawString(label, x + 20, y + dh);
        x += labelWidth + 60;

        label = "Retained";
        labelWidth = (int) fm.getStringBounds(label, g2D).getWidth();
        g2D.setColor(colorScale.getMidColor());
        g2D.fillRect(x, y, 10, 10);
        g2D.setColor(Color.BLACK);
        g2D.drawRect(x, y, 10, 10);
        g2D.drawString(label, x + 20, y + dh);
        x += labelWidth + 60;

        label = "Conflict";
        labelWidth = (int) fm.getStringBounds(label, g2D).getWidth();
        g2D.setColor(colorScale.getMinColor());
        g2D.fillRect(x, y, 10, 10);
        g2D.setColor(Color.BLACK);
        g2D.drawRect(x, y, 10, 10);
        g2D.drawString(label, x + 20, y + dh);
        x += labelWidth + 60;

    }
}
