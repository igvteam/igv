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

import org.broad.igv.feature.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
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

    private List<RegionFrame> renderedFrames = new ArrayList<RegionFrame>(1000);

    public CursorTrackPanel(CursorTrack t, CursorModel cursorModel, CursorMainPanel mainPanel) {
        this.track = t;
        this.model = cursorModel;
        this.mainPanel = mainPanel;
        MouseAdapter ma = new MouseHandler();
        addMouseListener(ma);
        addMouseMotionListener(ma);
        setToolTipText("");
    }

    public void setTrack(CursorTrack track) {
        this.track = track;
    }


    private void showPopup(int x, int y) {

        JPopupMenu menu = new JPopupMenu(); //track.getName());
//
//        JMenuItem setFramesItem = new JMenuItem("Check sort");
//        setFramesItem.addActionListener(new ActionListener() {
//            @Override
//            public void actionPerformed(ActionEvent e) {
//                checkSort();
//            }
//        });
//        menu.add(setFramesItem);
        CursorTrack.SignalField[] signalFields = track.getAllSignalFields();
        if (signalFields.length > 1) {
            ButtonGroup group = new ButtonGroup();
            for (CursorTrack.SignalField sf : signalFields) {
                final JRadioButton button = new JRadioButton(sf.toString());
                group.add(button);
                menu.add(button);
                if (sf == track.getSignalField()) button.setSelected(true);

                button.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        CursorTrack.SignalField sf = CursorTrack.SignalField.valueOf(button.getText());
                        track.setSignalField(sf);
                        mainPanel.repaint();
                    }
                });
            }
        }


        menu.show(this, x, y);

    }


    @Override
    protected void paintComponent(Graphics graphics) {

        super.paintComponent(graphics);

        renderedFrames.clear();

        if (model == null || track == null) return;

        List<CursorRegion> regionsList = model.getFilteredRegions();
        if (regionsList == null) return;

        int frameMargin = model.getFrameMargin();
        if (frameMargin < 1) {
            graphics.setColor(Color.white);
            graphics.fillRect(0, 0, getWidth(), getHeight());
        }

        Color c = track.getColor();
        graphics.setColor(c);

        double framePixelWidth = model.getFramePixelWidth();
        int frameBPWidth = model.getFrameBPWidth();
        double scale = ((double) frameBPWidth) / (framePixelWidth - frameMargin);
        double origin = model.getOrigin();   // In frame units
        double end = origin + ((double) getWidth()) / framePixelWidth;


        int startRegionNumber = (int) origin;
        if (startRegionNumber >= regionsList.size()) return;

        int h = getHeight();
        int sampleInterval = Math.max(1, (int) Math.round(1.0 / framePixelWidth));

        double max = track.getScale().getMax();

        for (int regionNumber = startRegionNumber; regionNumber < regionsList.size(); regionNumber += sampleInterval) {

            if (regionNumber > end) break;     // Absolutely critical for performance

            CursorRegion region = regionsList.get(regionNumber);
            String chr = region.getChr();

            int pxStart = (int) ((regionNumber - origin) * framePixelWidth + frameMargin / 2);
            int pxEnd = (int) (framePixelWidth > 1 ? pxStart + framePixelWidth :
                    (regionNumber + 1 - origin) * frameBPWidth) - frameMargin;

            renderedFrames.add(new RegionFrame(pxStart, pxEnd, region));

            int maxFeatureHeight = h - 10;

            // Region block
            if (frameMargin > 0) {
                int pxWidth = pxEnd - pxStart;
                graphics.setColor(Color.white);
                graphics.fillRect(pxStart, 0, pxWidth, h);
            }

            List<BasicFeature> features = track.getFeatures(chr);
            if (features == null) continue;

            int l2 = track.getLongestFeatureLength(chr);
            Iterator<BasicFeature> regionFeatures = region.getFeatureIterator(features, l2, model.getFrameBPWidth());

            double bpStart = region.getLocation() - frameBPWidth / 2;
            while (regionFeatures.hasNext()) {

                BasicFeature feature = regionFeatures.next();
                float score = track.getSignal(feature);
                if (Float.isNaN(score)) continue;

                int pStart = (int) Math.min(pxEnd, Math.max(pxStart, pxStart + (feature.getStart() - bpStart) / scale));
                int pEnd = (int) Math.min(pxEnd, pxStart + (feature.getEnd() - bpStart) / scale);
                int pw = Math.max(1, pEnd - pStart);

                graphics.setColor(c);

                // Height proportional to score
                int fh = (int) ((score / max) * maxFeatureHeight);
                graphics.fillRect(pStart, h - fh, pw, fh);

            }
        }


        graphics.setColor(Color.black);
        graphics.drawString(track.getName(), 10, 12);


    }

    @Override
    public String getToolTipText(MouseEvent event) {

        if (renderedFrames != null) {
            for (RegionFrame f : renderedFrames) {
                final int x = event.getX();
                if (x > f.pxStart && x <= f.pxEnd) {

                    int bpw = model.getFrameBPWidth();
                    int center = f.region.getLocation();
                    int start = center - bpw / 2;
                    int end = center + bpw / 2;
                    StringBuffer sb = new StringBuffer("<html>");
                    sb.append("Region " + f.region.getChr() + ":" + start + "-" + end);

                    String chr = f.region.chr;
                    List<BasicFeature> features = track.getFeatures(chr);
                    if (features != null) {
                        int longest = track.getLongestFeatureLength(chr);
                        Iterator<BasicFeature> iter = f.region.getFeatureIterator(features, longest, bpw);
                        while (iter.hasNext()) {
                            BasicFeature rf = iter.next();
                            sb.append("<br/>" + rf.getChr() + ":" + rf.getStart() + "-" + rf.getEnd() + " : " +
                                    track.getSignal(rf));
                        }
                    }

                    return sb.toString();
                }
            }
        }
        return "";
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


    static class RegionFrame {
        int pxStart;
        int pxEnd;
        CursorRegion region;

        RegionFrame(int pxStart, int pxEnd, CursorRegion region) {
            this.pxStart = pxStart;
            this.pxEnd = pxEnd;
            this.region = region;
        }
    }
}
