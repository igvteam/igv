package org.broad.igv.cursor;

import com.jidesoft.swing.JideButton;
import org.broad.igv.feature.*;
import org.broad.igv.ui.color.ColorUtilities;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.Serializable;
import java.util.List;

/**
 * @author jrobinso
 *         Date: 1/14/14
 *         Time: 12:41 PM
 */
public class CursorTrackPanel extends JComponent implements Serializable {

    Cursor model;
    CursorTrack track;
    private CursorMainPanel mainPanel;

    public CursorTrackPanel(CursorTrack t, Cursor cursor, CursorMainPanel mainPanel) {
        this.track = t;
        this.model = cursor;
        this.mainPanel = mainPanel;
        MouseAdapter ma = new MouseHandler();
        addMouseListener(ma);
        addMouseMotionListener(ma);
    }

    public static float getAlpha(float minRange, float maxRange, float value) {
        float binWidth = (maxRange - minRange) / 9;
        int binNumber = (int) ((value - minRange) / binWidth);
        return Math.min(1.0f, 0.2f + (binNumber * 0.8f) / 9);
    }


    public void setMainPanel(CursorMainPanel mainPanel) {
        this.mainPanel = mainPanel;
    }

    public CursorMainPanel getMainPanel() {
        return mainPanel;
    }

    public void setTrack(CursorTrack track) {
        this.track = track;
    }


    @Override
    protected void paintComponent(Graphics graphics) {

        super.paintComponent(graphics);

        if (model == null || track == null) return;

        List<CursorFrame> frameList = model.getFrames();
        if (frameList == null) return;

        int frameMargin = model.getFrameMargin();

        int framePixelWidth = model.getFramePixelWidth();
        int frameBPWidth = model.getFrameBPWidth();
        double scale = ((double) frameBPWidth) / (framePixelWidth - frameMargin);
        double origin = model.getOrigin();   // In frame units
        double end = origin + ((double) getWidth()) / framePixelWidth;

        int h = getHeight();


        int startFrameNumber = (int) origin;
        if (startFrameNumber >= frameList.size()) return;

        for (int frameNumber = startFrameNumber; frameNumber < frameList.size(); frameNumber++) {

            if (frameNumber > end) break;

            CursorFrame frame = frameList.get(frameNumber);
            String chr = frame.getChr();

            double bpStart = frame.getLocation() - frameBPWidth / 2;
            double bpEnd = frame.getLocation() + frameBPWidth / 2;
            double pxStart = (frameNumber - origin) * framePixelWidth;
            double pxEnd = pxStart + framePixelWidth;

            int maxFeatureHeight = h;

            // Region block

            graphics.setColor(Color.white);
            graphics.fillRect((int) (pxStart + frameMargin / 2), 0, framePixelWidth - frameMargin, h);

            List<BasicFeature> features = track.getFeatures(chr);
            if (features == null) continue;

            int i0 = FeatureUtils.getIndexBefore(bpStart, features);
            if (i0 < 0) continue;

            for (int fIdx = i0; fIdx < features.size(); fIdx++) {

                LocusScore f = features.get(fIdx);
                if (f.getStart() > bpEnd) break;
                else if (f.getEnd() < bpStart) continue;
                else {

                    int pStart = (int) Math.max(pxStart + frameMargin / 2, pxStart + frameMargin / 2 + (f.getStart() - bpStart) / scale);
                    int pEnd = (int) Math.min(pxEnd - frameMargin / 2, pxStart + frameMargin / 2 + (f.getEnd() - bpStart) / scale);

                    Color c = track.getColor();
                    float min = 0;
                    float max = 1000;

                    float score = f.getScore();
                    float alpha = Float.isNaN(score) ? 1 : getAlpha(min, max, score);
                    c = ColorUtilities.getCompositeColor(c, alpha);
                    graphics.setColor(c);

                    int pw = pEnd == pStart ? 1 : pEnd - pStart;

                    // Height proportional to score
                    int fh = Math.max(10, (int) ((f.getScore() / 1000) * maxFeatureHeight));
                    graphics.fillRect(pStart, h - fh, pw, fh);

                }
            }

        }
    }


    class MouseHandler extends MouseAdapter {

        int mouseX;

        @Override
        public void mousePressed(MouseEvent mouseEvent) {
            mouseX = mouseEvent.getX();
        }

        @Override
        public void mouseDragged(MouseEvent mouseEvent) {
            int delta = mouseX - mouseEvent.getX();
            model.shiftOriginPixels(delta);
            mouseX = mouseEvent.getX();

            // TODO -- use event here and remove back pointer

            mainPanel.repaint();

        }
    }


}
