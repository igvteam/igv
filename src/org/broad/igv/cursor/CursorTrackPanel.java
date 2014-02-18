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

    private List<RegionFrame> renderedFrames;

    public CursorTrackPanel(CursorTrack t, CursorModel cursorModel, CursorMainPanel mainPanel) {
        this.track = t;
        this.model = cursorModel;
        this.mainPanel = mainPanel;
        MouseAdapter ma = new MouseHandler();
        addMouseListener(ma);
        addMouseMotionListener(ma);
        setToolTipText("");
        ToolTipManager.sharedInstance().setDismissDelay(20000);
    }

    public void setTrack(CursorTrack track) {
        this.track = track;
    }


    private void showPopup(int x, int y) {

        JPopupMenu menu = new JPopupMenu(); //track.getName());

        JMenuItem setFramesItem = new JMenuItem("Check sort");
        setFramesItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                checkSort();
            }
        });
        menu.add(setFramesItem);

        menu.show(this, x, y);

    }


    @Override
    protected void paintComponent(Graphics graphics) {

        super.paintComponent(graphics);

        renderedFrames = new ArrayList<RegionFrame>(1000);

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
        for (int regionNumber = startRegionNumber; regionNumber < regionsList.size(); regionNumber += sampleInterval) {

            if (regionNumber > end) break;     // Absolutely critical for performance

            CursorRegion region = regionsList.get(regionNumber);
            String chr = region.getChr();

            double bpStart = region.getLocation() - frameBPWidth / 2;
            double bpEnd = region.getLocation() + frameBPWidth / 2;

            int pxStart = (int) ((regionNumber - origin) * framePixelWidth);
            int pxEnd = (int) (pxStart + framePixelWidth);
            renderedFrames.add(new RegionFrame(pxStart, pxEnd, region));

            int maxFeatureHeight = h - 10;
            float max = 1000;

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

            while (regionFeatures.hasNext()) {
                BasicFeature feature = regionFeatures.next();
                float score = feature.getScore();
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
                if (event.getX() > f.pxStart && event.getX() <= f.pxEnd) {

                    int bpw = model.getFrameBPWidth();
                    int center = f.region.getLocation();
                    int start = center - bpw / 2;
                    int end = center + bpw / 2;
                    StringBuffer sb = new StringBuffer("<html>");
                    sb.append("Region " + f.region.getChr() + ":" + start + "-" + end + " : " + f.region.getTmp());

                    String chr = f.region.chr;
                    List<BasicFeature> features = track.getFeatures(chr);
                    if (features != null) {
                        int longest = track.getLongestFeatureLength(chr);
                        Iterator<BasicFeature> iter = f.region.getFeatureIterator(features, longest, bpw);
                        while (iter.hasNext()) {
                            BasicFeature rf = iter.next();
                            sb.append("<br/>" + rf.getChr() + ":" + rf.getStart() + "-" + rf.getEnd() + " : " + rf.getScore());
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


    protected void checkSort() {


        if (model == null || track == null) return;

        List<CursorRegion> regionsList = model.getFilteredRegions();
        if (regionsList == null) return;


        double framePixelWidth = model.getFramePixelWidth();
        int frameBPWidth = model.getFrameBPWidth();
        double origin = model.getOrigin();   // In frame units
        double end = origin + ((double) getWidth()) / framePixelWidth;


        int startRegionNumber = (int) origin;
        if (startRegionNumber >= regionsList.size()) return;

        int sampleInterval = Math.max(1, (int) Math.round(1.0 / framePixelWidth));
        double lastScore = -Double.MAX_VALUE;
        for (int regionNumber = startRegionNumber; regionNumber < regionsList.size(); regionNumber += sampleInterval) {

            if (regionNumber > end) break;     // Absolutely critical for performance

            CursorRegion region = regionsList.get(regionNumber);
            String chr = region.getChr();
            // Region block

            List<BasicFeature> features = track.getFeatures(chr);
            int l1 = track.getLongestFeatureLength(region.getChr());
            double score = region.getScore(features, l1, frameBPWidth);
            if (score < lastScore) {
                System.out.println(score);
            }
            lastScore = score;
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
