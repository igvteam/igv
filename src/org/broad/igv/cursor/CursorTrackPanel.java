package org.broad.igv.cursor;

import org.broad.igv.feature.*;
import org.broad.igv.ui.color.ColorUtilities;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
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

    CursorModel model;
    CursorTrack track;
    private CursorMainPanel mainPanel;

    public CursorTrackPanel(CursorTrack t, CursorModel cursorModel, CursorMainPanel mainPanel) {
        this.track = t;
        this.model = cursorModel;
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


    private void showPopup(int x, int y) {

        JPopupMenu menu = new JPopupMenu(); //track.getName());

        JMenuItem setFramesItem = new JMenuItem("Set frames");
        setFramesItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                model.setFrames(CursorUtils.createRegions(track));
                mainPanel.repaint();
            }
        });
        menu.add(setFramesItem);

        menu.show(this, x, y);

    }


    @Override
    protected void paintComponent(Graphics graphics) {

        super.paintComponent(graphics);

        if (model == null || track == null) return;

        List<CursorRegion> frameList = model.getFrames();
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

            CursorRegion frame = frameList.get(frameNumber);
            String chr = frame.getChr();

            double bpStart = frame.getLocation() - frameBPWidth / 2;
            double bpEnd = frame.getLocation() + frameBPWidth / 2;
            double pxStart = (frameNumber - origin) * framePixelWidth;
            double pxEnd = pxStart + framePixelWidth;

            int maxFeatureHeight = h-10;

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
                    if (alpha < 1) {
                        c = ColorUtilities.getCompositeColor(c, alpha);
                    }
                    graphics.setColor(c);

                    int pw = pEnd == pStart ? 1 : pEnd - pStart;

                    // Height proportional to score
                    int fh = Math.max(10, (int) ((f.getScore() / 1000) * maxFeatureHeight));
                    graphics.fillRect(pStart, h - fh, pw, fh);

                }
            }

            graphics.setColor(Color.black);
            graphics.drawString(track.getName(), 10, 12);

        }
    }


    class MouseHandler extends MouseAdapter {

        int mouseX;

        @Override
        public void mousePressed(MouseEvent e) {
            mouseX = e.getX();
            evaluatePopup(e);
        }

        @Override
        public void mouseReleased(MouseEvent e) {
            evaluatePopup(e);
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


    private void evaluatePopup(MouseEvent e) {
        if (e.isPopupTrigger()) {
            showPopup(e.getX(), e.getY());
        }
    }


}
