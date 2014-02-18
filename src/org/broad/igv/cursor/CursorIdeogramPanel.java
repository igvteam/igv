package org.broad.igv.cursor;

import org.broad.igv.feature.BasicFeature;
import org.broad.igv.ui.color.ColorUtilities;

import javax.swing.*;
import java.awt.*;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * @author jrobinso
 *         Date: 1/15/14
 *         Time: 9:22 AM
 */
public class CursorIdeogramPanel extends JComponent implements Serializable {

    CursorModel model;
    List<CursorTrack> tracks = new ArrayList<CursorTrack>();
    boolean drawViewRect = true;

    public CursorIdeogramPanel() {
        //     setBorder(BorderFactory.createLineBorder(Color.black));

    }

    public static float getAlpha(float minRange, float maxRange, float value) {
        float binWidth = (maxRange - minRange) / 9;
        int binNumber = (int) ((value - minRange) / binWidth);
        return Math.min(1.0f, 0.2f + (binNumber * 0.8f) / 9);
    }

    @Override
    protected void paintComponent(Graphics graphics) {

        super.paintComponent(graphics);

        graphics.setColor(Color.gray);
        graphics.drawLine(0, 0, 0, getHeight());

        if (tracks.size() > 0) {
            paintBackground(graphics);
        }

        if (model != null && model.getFilteredRegions() != null && drawViewRect) {

            int length = model.getFilteredRegions().size();
            if (length == 0) return;

            int px = (int) ((model.getOrigin() / length) * getWidth());

            double nFrames = ((double) getWidth()) / model.getFramePixelWidth();
            int width = Math.max(1, (int) ((nFrames / length) * getWidth()));

            graphics.setColor(Color.black);
            graphics.drawRect(px, 0, width, getHeight() - 1);
            graphics.drawRect(px + 1, 1, width - 2, getHeight() - 2);
        }

    }

    public void setModel(CursorModel model) {
        this.model = model;
    }


    private void paintBackground(Graphics graphics) {

        if (model == null || tracks.isEmpty()) return;

        List<CursorRegion> frameList = model.getFilteredRegions();
        if (frameList == null) return;

        // We'll sample frames and give each 1 pixel
        double sampleInterval = ((double) frameList.size()) / getWidth();

        int bh = getHeight() - 2;
        double dh = ((double) bh) / tracks.size();
        int px = 0;
        for (double frameNumber = 0; frameNumber < frameList.size(); frameNumber += sampleInterval) {

            CursorRegion frame = frameList.get((int) frameNumber);
            String chr = frame.getChr();
            int maxFeatureHeight = (int) dh;

            graphics.setColor(Color.white);
            graphics.drawLine(px, 0, px, getHeight());

            double base = dh + 1;
            for (CursorTrack track : tracks) {

                List<BasicFeature> features = track.getFeatures(chr);
                if (features == null) continue;

                int l2 = track.getLongestFeatureLength(chr);
                Iterator<BasicFeature> regionFeatures = frame.getFeatureIterator(features, l2, model.getFrameBPWidth());

                while (regionFeatures.hasNext()) {
                    BasicFeature f = regionFeatures.next();

                    Color c = track.getColor();
                    float min = 0;
                    float max = 1000;

                    float score = f.getScore();
                    float alpha = Float.isNaN(score) ? 1 : getAlpha(min, max, score);
                    c = ColorUtilities.getCompositeColor(c, alpha);
                    graphics.setColor(c);

                    graphics.drawLine(px, (int) base - maxFeatureHeight, px, (int) base);

                }
                base += dh;
            }
            px++;
        }
    }

    public void addTrack(CursorTrack track) {
        this.tracks.add(track);
    }

}
