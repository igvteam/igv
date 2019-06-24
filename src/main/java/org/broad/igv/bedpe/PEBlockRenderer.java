package org.broad.igv.bedpe;

import org.broad.igv.track.RenderContext;

import java.awt.*;
import java.util.List;

public class PEBlockRenderer implements BedPERenderer {

    BedPETrack track;
    int rowHeight = 10;

    public PEBlockRenderer(BedPETrack track) {
        this.track = track;
    }

    @Override
    public void render(List<BedPE> features, RenderContext context, Rectangle trackRectangle) {

        Graphics2D g = null;

        try {
            g = (Graphics2D) context.getGraphics().create();
            String chr = context.getChr();
            double origin = context.getOrigin();
            double locScale = context.getScale();
            Color trackColor = track.getColor();

            final int baseY = trackRectangle.y;

            for (BedPE bedPE : features) {

                BedPEFeature feature = bedPE.get();
                Color fcolor = feature.color == null ? trackColor : feature.color;
                if (fcolor != null) {
                    g.setColor(fcolor);
                }

                int blockY = baseY + bedPE.getRow() * rowHeight;

                if (feature.isSameChr()) {
                    int ps1 = (int) ((feature.start1 - origin) / locScale);
                    int pe1 = (int) ((feature.end1 - origin) / locScale);
                    if (pe1 >= trackRectangle.getX() && ps1 <= trackRectangle.getMaxX()) {
                        ps1 = drawBlock(ps1, pe1, baseY, g);
                    }

                    int ps2 = (int) ((feature.start2 - origin) / locScale);
                    int pe2 = (int) ((feature.end2 - origin) / locScale);
                    if (pe2 >= trackRectangle.getX() && ps2 <= trackRectangle.getMaxX()) {
                        ps2 = drawBlock(ps2, pe2, baseY, g);
                    }

                    // connecting line
                    if (feature.isSameChr()) {
                        int pl1 = Math.min(pe1, pe2);
                        int pl2 = Math.max(ps1, ps2);
                        final int connectorY = baseY + rowHeight / 2;
                        g.drawLine(pl1, connectorY, pl2, connectorY);
                    }
                } else {
                    int ps1 = (int) ((feature.getStart() - origin) / locScale);
                    int pe1 = (int) ((feature.getEnd() - origin) / locScale);
                    if (pe1 >= trackRectangle.getX() && ps1 <= trackRectangle.getMaxX()) {
                        drawBlock(ps1, pe1, baseY, g);
                    }

                }
            }
        } finally {
            if (g != null) g.dispose();
        }
    }

    private int drawBlock(int ps1, int pe1, int blockY, Graphics2D g) {
        // Trim width if possible to insure a gap between blocks
        int w1 = Math.max(1, pe1 - ps1);
        if (w1 > 3) w1--;
        if (w1 > 5) ps1++;
        g.fillRect(ps1, blockY, w1, rowHeight);
        return ps1;
    }

}
