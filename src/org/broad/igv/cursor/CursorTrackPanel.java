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
                model.setRegions(CursorUtils.createRegions(track));
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

        List<CursorRegion> regionsList = model.getFilteredRegions();
        if (regionsList == null) return;

        int frameMargin = model.getFrameMargin();

        double framePixelWidth = model.getFramePixelWidth();
        int frameBPWidth = model.getFrameBPWidth();
        double scale = ((double) frameBPWidth) / (framePixelWidth - frameMargin);
        double origin = model.getOrigin();   // In frame units
        double end = origin + ((double) getWidth()) / framePixelWidth;


        int startRegionNumber = (int) origin;
        if (startRegionNumber >= regionsList.size()) return;

        int h = getHeight();
        int sampleInterval = Math.max(1, (int) Math.round(1.0 / framePixelWidth));
        for (int regionNumber = startRegionNumber; regionNumber < regionsList.size(); regionNumber += sampleInterval) {

            if (regionNumber > end) break;     // Absolutely critical for performance

            CursorRegion region = regionsList.get(regionNumber);
            String chr = region.getChr();

            double bpStart = region.getLocation() - frameBPWidth / 2;
            double bpEnd = region.getLocation() + frameBPWidth / 2;
            int pxStart = (int) ((regionNumber - origin) * framePixelWidth + frameMargin / 2);
            int pxEnd = (int) (framePixelWidth > 1 ? pxStart + framePixelWidth :
                    (regionNumber + 1 - origin) * frameBPWidth) - frameMargin;
            int pxWidth = pxEnd - pxStart;

            int maxFeatureHeight = h - 10;

            // Region block

            graphics.setColor(Color.white);
            graphics.fillRect(pxStart, 0, pxWidth, h);

            List<BasicFeature> features = track.getFeatures(chr);
            if (features == null) continue;


            int l2 = track.getLongestFeatureLength(chr);
            double s0 = l2 < 0 ? 0 : bpStart - l2;
            int i0 = FeatureUtils.getIndexBefore(s0, features);
            if (i0 < 0) continue;

            for (int fIdx = i0; fIdx < features.size(); fIdx++) {

                LocusScore feature = features.get(fIdx);
                if (feature.getStart() >= bpEnd) break;
                else if (feature.getEnd() > bpStart && feature.getStart() <= bpEnd) {


                    /*// For debugging.
                    int l1 = track.getLongestFeatureLength(region.getChr());
                    float score = (float) region.getScore(features, l1, frameBPWidth);
                    int pStart = pxStart;
                    int pEnd = pxEnd;*/

                    float score = feature.getScore();
                    int pStart = (int) Math.max(pxStart, pxStart + (feature.getStart() - bpStart) / scale);
                    int pEnd = (int) Math.min(pxEnd, pxStart + (feature.getEnd() - bpStart) / scale);

                    Color c = track.getColor();
                    float min = 0;
                    float max = 1000;

                   // float alpha = Float.isNaN(score) ? 1 : getAlpha(min, max, score);
                   // if (alpha < 1) {
                   //     c = ColorUtilities.getCompositeColor(c, alpha);
                   // }
                    graphics.setColor(c);

                    int pw = Math.max(1, pEnd - pStart);

                    // Height proportional to score
                    int fh = (int) ((score / max) * maxFeatureHeight);
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
