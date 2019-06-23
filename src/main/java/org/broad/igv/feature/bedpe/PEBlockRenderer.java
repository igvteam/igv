package org.broad.igv.feature.bedpe;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.RenderContext;

import java.awt.*;
import java.util.List;

public class PEBlockRenderer implements BedPERenderer {

    BedPETrack track;

    public PEBlockRenderer(BedPETrack track) {
        this.track = track;
    }

    @Override
    public void render(List<BedPE> features, RenderContext context, Rectangle trackRectangle) {

        Graphics2D g = null;

        try {
//            g = (Graphics2D) context.getGraphics().create();
//            String chr = context.getChr();
//            double origin = context.getOrigin();
//            double locScale = context.getScale();
//            Color trackColor = track.getColor();
//
//            for (BedPE bedPE : features) {
//
//                BedPEFeature feature = bedPE.get();
//                Color fcolor = feature.color == null ? trackColor : feature.color;
//                if (fcolor != null) {
//                    g.setColor(fcolor);
//                }
//
//                final int h = 10;
//                final int blockY = trackRectangle.y + trackRectangle.height - h;
//
//                int ps1 = (int) ((feature.start1 - origin) / locScale);
//                int pe1 = (int) ((feature.end1 - origin) / locScale);
//
//                if (feature.isSameChr()) {
//                    String chr1 = genome == null ? feature.chr1 : genome.getCanonicalChrName(feature.chr1);
//                    if (chr1.equals(chr)) {
//                        // Trim width if possible to insure a gap between blocks
//                        int w1 = Math.max(1, pe1 - ps1);
//                        if (w1 > 3) w1--;
//                        if (w1 > 5) ps1++;
//                        if (pe1 >= trackRectangle.getX() && ps1 <= trackRectangle.getMaxX()) {
//                            g.fillRect(ps1, blockY, w1, 10);
//                        }
//                    }
//
//
//                    int ps2 = (int) ((feature.start2 - origin) / locScale);
//                    int pe2 = (int) ((feature.end2 - origin) / locScale);
//                    String chr2 = genome == null ? feature.chr2 : genome.getCanonicalChrName(feature.chr2);
//                    if (chr2.equals(chr)) {
//                        // Trim width if possible to insure a gap between blocks
//                        int w2 = Math.max(1, pe2 - ps2);
//                        if (w2 > 3) w2--;
//                        if (w2 > 5) ps2++;
//
//                        if (pe2 >= trackRectangle.getX() && ps2 <= trackRectangle.getMaxX()) {
//                            g.fillRect(ps2, blockY, w2, 10);
//                        }
//                    }
//
//                    // connecting line
//                    if (feature.isSameChr()) {
//                        int pl1 = Math.min(pe1, pe2);
//                        int pl2 = Math.max(ps1, ps2);
//                        final int connectorY = blockY + h / 2;
//                        g.drawLine(pl1, connectorY, pl2, connectorY);
//                    }
//                }
//            }
        } finally {
            if (g != null) g.dispose();
        }
    }

}
