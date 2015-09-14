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

package org.broad.igv.peaks;

import org.apache.log4j.Logger;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackClickEvent;
import org.broad.igv.track.TrackMenuUtils;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.lang.ref.SoftReference;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * @author jrobinso
 * @date Apr 23, 2011
 */
public class PeakTrackMenu extends IGVPopupMenu {

    private static Logger log = Logger.getLogger(PeakTrackMenu.class);
    private PeakTrack track;
    private JRadioButtonMenuItem colorByScoreMI;
    private JRadioButtonMenuItem colorByFoldMI;

    public PeakTrackMenu(PeakTrack track, TrackClickEvent te) {
        super();
        // TODO -- what if multiple tracks are selected?
        this.track = track;
        init(te);
    }

    private void init(final TrackClickEvent trackClickEvent) {

        Collection<Track> tracks = IGV.getInstance().getSelectedTracks(); //Arrays.asList(new Track[] {track});

        //Title
        JLabel popupTitle = new JLabel("<html><b>" + track.getName(), JLabel.LEFT);
        Font newFont = getFont().deriveFont(Font.BOLD, 12);
        popupTitle.setFont(newFont);
        add(popupTitle);
        addSeparator();

        add(TrackMenuUtils.getTrackRenameItem(tracks));
        add(TrackMenuUtils.getChangeTrackHeightItem(tracks));
        add(TrackMenuUtils.getChangeFontSizeItem(tracks));

        //Change Track Settings
        addDisplayModeItems();

        addSeparator();
        JMenuItem plotItem = new JMenuItem("Open Trend Plot...");
        plotItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent actionEvent) {
                openTimeSeriesPlot(trackClickEvent);
            }
        });
        add(plotItem);

        if (track.isShowSignals()) {
            add(TrackMenuUtils.getDataRangeItem(tracks));
            add(TrackMenuUtils.getLogScaleItem(tracks));
            add(TrackMenuUtils.getShowDataRangeItem(tracks));
        }

        addSeparator();
        add(TrackMenuUtils.getRemoveMenuItem(tracks));
    }

    public void addDisplayModeItems() {
        addSeparator();
        ButtonGroup group = new ButtonGroup();

        Track.DisplayMode displayMode = track.getDisplayMode();

        JRadioButtonMenuItem m1 = new JRadioButtonMenuItem("Compressed");
        m1.setSelected(displayMode == Track.DisplayMode.COLLAPSED);
        m1.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                track.setDisplayMode(Track.DisplayMode.COLLAPSED);
                PeakTrack.setAnimate(false);
                IGV.getInstance().doRefresh();
            }
        });

        JRadioButtonMenuItem m3 = new JRadioButtonMenuItem("Time Series");
        m3.setSelected(displayMode == Track.DisplayMode.EXPANDED);
        m3.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                track.setDisplayMode(Track.DisplayMode.EXPANDED);
                PeakTrack.setAnimate(false);
                IGV.getInstance().doRefresh();
            }
        });

        JRadioButtonMenuItem m4 = new JRadioButtonMenuItem("Animate");
        m4.setSelected(PeakTrack.animate);
        m4.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                track.setDisplayMode(Track.DisplayMode.COLLAPSED);
                PeakTrack.setAnimate(true);
             }
        });

        add(m1);
        add(m3);
        add(m4);
        group.add(m1);
        group.add(m3);
        group.add(m4);

    }

    public void addColorByItems() {


        addSeparator();
        add(new JLabel("Color by"));

        ButtonGroup group = new ButtonGroup();

        colorByScoreMI = new JRadioButtonMenuItem("Score");
        colorByScoreMI.setSelected(PeakTrack.getColorOption() == PeakTrack.ColorOption.SCORE);
        colorByScoreMI.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                PeakTrack.setShadeOption(PeakTrack.ColorOption.SCORE);
                IGV.getInstance().repaint();
            }
        });

        colorByFoldMI = new JRadioButtonMenuItem("Fold change");
        colorByFoldMI.setSelected(PeakTrack.getColorOption() == PeakTrack.ColorOption.FOLD_CHANGE);
        colorByFoldMI.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                PeakTrack.setShadeOption(PeakTrack.ColorOption.FOLD_CHANGE);
                IGV.getInstance().repaint();
            }
        });

        add(colorByScoreMI);
        add(colorByFoldMI);
        group.add(colorByScoreMI);
        group.add(colorByFoldMI);

    }


    void openTimeSeriesPlot(TrackClickEvent trackClickEvent) {

        if (trackClickEvent == null) return;

        ReferenceFrame referenceFrame = trackClickEvent.getFrame();
        if (referenceFrame == null) return;

        String chr = referenceFrame.getChrName();
        double position = trackClickEvent.getChromosomePosition();

        XYSeriesCollection data = new XYSeriesCollection();
        List<Color> colors = new ArrayList();
        for (SoftReference<PeakTrack> ref : PeakTrack.instances) {
            PeakTrack track = ref.get();
            if (track != null) {
                Peak peak = track.getFilteredPeakNearest(chr, position);
                if (peak != null) {
                    XYSeries series = new XYSeries(track.getName());
                    float[] scores = peak.getTimeScores();
                    if (scores.length == 4) {
                        float t0 = scores[0] + 10;

                        series.add(0, (scores[0] + 10) / t0);
                        series.add(30, (scores[1] + 10) / t0);
                        series.add(60, (scores[2] + 10) / t0);
                        series.add(120, (scores[3] + 10) / t0);
                    }
                    data.addSeries(series);
                    Color c = track.getName().contains("Pol") ? Color.black : track.getColor();
                    colors.add(c);
                }
            }
        }

        final JFreeChart chart = ChartFactory.createXYLineChart(
                "",
                "Time",
                "Score",
                data,
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );

        NumberAxis axis = (NumberAxis) chart.getXYPlot().getDomainAxis(0);
        axis.setTickUnit(new NumberTickUnit(30));

        final ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(400, 400));
        chartPanel.setSize(new java.awt.Dimension(400, 400));


        XYItemRenderer renderer = chart.getXYPlot().getRenderer();
        for (int i = 0; i < colors.size(); i++) {
            renderer.setSeriesPaint(i, colors.get(i));
        }

        chartPanel.setBackground(Color.white);
        chart.getXYPlot().setBackgroundPaint(Color.white);
        chart.getXYPlot().setRangeGridlinePaint(Color.black);


        PeakTimePlotFrame frame = new PeakTimePlotFrame(chartPanel);
        frame.setVisible(true);


    }
}
