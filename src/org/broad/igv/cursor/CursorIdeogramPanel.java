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

    public static double getAlpha(double minRange, double maxRange, double value) {
        double binWidth = (maxRange - minRange) / 9;
        int binNumber = (int) ((value - minRange) / binWidth);
        return Math.min(1.0, 0.2 + (binNumber * 0.8) / 9);
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

                CursorTrack.Range yScale = track.getScale();
                double min = yScale.getMin();
                double max = yScale.getMax();

                List<BasicFeature> features = track.getFeatures(chr);
                if (features == null) continue;

                int l2 = track.getLongestFeatureLength(chr);
                Iterator<BasicFeature> regionFeatures = frame.getFeatureIterator(features, l2, model.getFrameBPWidth());

                while (regionFeatures.hasNext()) {
                    BasicFeature f = regionFeatures.next();

                    Color c = track.getColor();

                    float score = track.getSignal(f);
                    double alpha = Float.isNaN(score) ? 1 : getAlpha(min, max, score);
                    c = ColorUtilities.getCompositeColor(c, (float) alpha);
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
