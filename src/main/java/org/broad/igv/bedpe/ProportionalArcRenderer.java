package org.broad.igv.bedpe;

import org.broad.igv.track.RenderContext;

import java.awt.*;
import java.awt.geom.Arc2D;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.broad.igv.bedpe.BedPETrack.Direction.UP;

public class ProportionalArcRenderer implements BedPERenderer {


    private Map<Color, Color> alphaColors = new HashMap<>();

    BedPETrack track;
    private double logMaxScore = 1;

    public ProportionalArcRenderer(BedPETrack track) {
        this.track = track;
    }

    public void render(List<BedPE> features, RenderContext context, Rectangle trackRectangle) {

        Graphics2D g = null;

        try {
            g = (Graphics2D) context.getGraphics().create();
            double origin = context.getOrigin();
            double locScale = context.getScale();
            Color trackColor = track.getColor();

            if (track.thickness > 1) {
                g.setStroke(new BasicStroke(track.thickness));
            }

            if (track.autoscale) {
                autoscale(features);
            }

            for (BedPE bedPE : features) {

                double p1 = (bedPE.getStart() - origin) / locScale;
                double p2 = (bedPE.getEnd() - origin) / locScale;

                if (p2 >= trackRectangle.getX() && p1 <= trackRectangle.getMaxX()) {

                    BedPETrack.Direction direction = track.direction;
                    int gap = track.gap;
                    int h = trackRectangle.height - gap;
                    double logMax = logMaxScore;
                    if (logMax > 0 && bedPE.getScore() > 0) {
                        h = (int) ((Math.log10(bedPE.getScore()) / logMax) * h);
                    }

                    if (bedPE.isSameChr()) {

                        BedPEFeature feature = bedPE.get();
                        Color fcolor = feature.color == null ? trackColor : feature.color;
                        if (fcolor != null) {
                            g.setColor(fcolor);
                        }

                        double pixelStart = (feature.getMidStart() - origin) / locScale;
                        double pixelEnd = (feature.getMidEnd() - origin) / locScale;
                        int w = (int) (pixelEnd - pixelStart);
                        if (w < 3) {
                            w = 3;
                            pixelStart--;
                        }



                        double y = direction == UP ? gap + trackRectangle.y + trackRectangle.height - h : gap + trackRectangle.y - h;
                        int angleSt = direction == UP ? 0 : 180;
                        Arc2D.Double arcPath = new Arc2D.Double(
                                pixelStart,
                                y,
                                w,
                                2 * h,
                                angleSt,
                                180,
                                Arc2D.OPEN
                        );

                        g.draw(arcPath);
                        Color shadedColor = getAlphaColor(fcolor, 0.05f);
                        g.setColor(shadedColor);
                        g.fill(arcPath);

                    } else {
                        Color fcolor = bedPE.get().color == null ? Color.black : bedPE.get().color;
                        g.setColor(fcolor);
                        double ps = ((bedPE.getStart() + bedPE.getEnd()) / 2 - origin) / locScale;
                        int yBase = direction == UP ? trackRectangle.y + trackRectangle.height - h : trackRectangle.y + gap;
                        g.drawLine((int) ps, yBase, (int) ps, yBase + h);
                    }

                }
            }
        } finally {
            if (g != null) g.dispose();
        }
    }

    /**
     * Autoscale max height -- specific to proportional arc mode
     *
     * @param features
     */
    private void autoscale(List<BedPE> features) {
        double maxScore = 0;
        for (BedPE f : features) {
            maxScore = Math.max(maxScore, f.getScore());
        }
        if (maxScore > 0) {
            logMaxScore = Math.log10(maxScore);
        } else {
            logMaxScore = 1;
        }
    }

    private Color getAlphaColor(Color fcolor, float alpha) {
        Color ac = alphaColors.get(fcolor);
        if (ac == null) {
            float[] rgb = new float[3];
            rgb = fcolor.getRGBColorComponents(rgb);
            ac = new Color(rgb[0], rgb[1], rgb[2], alpha);
            alphaColors.put(fcolor, ac);
        }
        return ac;
    }

}

