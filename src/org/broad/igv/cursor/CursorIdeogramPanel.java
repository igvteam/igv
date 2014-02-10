package org.broad.igv.cursor;

import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.ui.color.ColorUtilities;

import javax.swing.*;
import java.awt.*;
import java.io.Serializable;
import java.util.ArrayList;
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
        setBorder(BorderFactory.createLineBorder(Color.black));
    }

    @Override
    protected void paintComponent(Graphics graphics) {

        super.paintComponent(graphics);

        if (tracks.size() > 0) {
            paintBackground(graphics);
        }

        if (model != null && model.getFilteredRegions() != null && drawViewRect) {

            int length = model.getFilteredRegions().size();
            if (length == 0) return;

            int px = (int) ((model.getOrigin() / length) * getWidth());

            double nFrames = ((double) getWidth()) / model.getFramePixelWidth();
            int width = Math.min(1, (int) ((nFrames / length) * getWidth()));

            graphics.drawRect(px, 1, width, getHeight() - 2);
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

        // TODO -- sampleInterval < 1;

        int frameBPWidth = model.getFrameBPWidth();
        double dh = ((double) getHeight()) / tracks.size();
        int px = 0;
        for (double frameNumber = 0; frameNumber < frameList.size(); frameNumber += sampleInterval) {

            CursorRegion frame = frameList.get((int) frameNumber);
            String chr = frame.getChr();

            int maxFeatureHeight = (int) dh;

            graphics.setColor(Color.white);
            graphics.drawLine(px, 0, px, getHeight());

            double base = dh;
            for (CursorTrack track : tracks) {

                List<BasicFeature> features = track.getFeatures(chr);
                if (features == null) continue;

                double bpStart = frame.getLocation() - frameBPWidth / 2;
                double bpEnd = frame.getLocation() + frameBPWidth / 2;
                int i0 = FeatureUtils.getIndexBefore(bpStart, features);
                if (i0 < 0) continue;

                for (int fIdx = i0; fIdx < features.size(); fIdx++) {

                    LocusScore f = features.get(fIdx);
                    if (f.getStart() > bpEnd) break;
                    else if (f.getEnd() < bpStart) continue;
                    else {


                        Color c = track.getColor();
                        float min = 0;
                        float max = 1000;

                        float score = f.getScore();
                        float alpha = Float.isNaN(score) ? 1 : CursorTrackPanel.getAlpha(min, max, score);
                        c = ColorUtilities.getCompositeColor(c, alpha);
                        graphics.setColor(c);


                        // Height proportional to score
                        int fh = f.getScore() == 0 ? 0 : Math.max(1, (int) ((f.getScore() / 1000) * maxFeatureHeight));
                        graphics.drawLine(px, (int) base - maxFeatureHeight, px, (int) base);

                    }
                }
                base += dh;
            }
            px++;
        }
    }

    public void addTrack(CursorTrack track) {
        this.tracks.add(track);
    }

    public void setDrawViewRect(boolean drawViewRect) {
        this.drawViewRect = drawViewRect;
    }
}
