package org.broad.igv.ui.panel;

import org.broad.igv.charts.ScatterPlotUtils;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.IGV;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Set;

/**
 * @author jrobinso
 * @date Sep 18, 2010
 */
public class RegionMenu extends JPopupMenu {

    RegionOfInterest roi;
    ReferenceFrame frame;

    public RegionMenu(final RegionOfInterest roi, final ReferenceFrame frame) {
        this(roi, frame, null);
    }

    public RegionMenu(final RegionOfInterest roi, final ReferenceFrame frame, String title) {

        this.roi = roi;
        this.frame = frame;

        if (title != null) {
            JLabel heading = new JLabel("<html>&nbsp;<b>" + title);
            //heading.setFont(UIConstants.boldFont);
            add(heading);
            addSeparator();
        }


        Set<TrackType> loadedTypes = IGV.getInstance().getLoadedTypes();

        if (loadedTypes.contains(TrackType.COPY_NUMBER) ||
                loadedTypes.contains(TrackType.ALLELE_SPECIFIC_COPY_NUMBER) ||
                loadedTypes.contains(TrackType.CNV)) {
            JMenuItem item = new JMenuItem("Sort by amplification");
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    IGV.getInstance().sortByRegionScore(roi, RegionScoreType.AMPLIFICATION, frame);
                }
            });
            add(item);


            item = new JMenuItem("Sort by deletion");
            item.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {

                    IGV.getInstance().sortByRegionScore(roi, RegionScoreType.DELETION, frame);
                    IGV.getInstance().getContentPane().repaint();
                }
            });
            add(item);

            item = new JMenuItem("Sort by breakpoint amplitudes");
            item.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {

                    IGV.getInstance().sortByRegionScore(roi, RegionScoreType.FLUX, frame);
                    IGV.getInstance().getContentPane().repaint();
                }
            });
            add(item);
        }

        if (loadedTypes.contains(TrackType.GENE_EXPRESSION)) {
            JMenuItem item = new JMenuItem("Sort by expression");
            item.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {

                    IGV.getInstance().sortByRegionScore(roi, RegionScoreType.EXPRESSION, frame);
                    IGV.getInstance().getContentPane().repaint();

                }
            });

            add(item);
        }

        if (loadedTypes.contains(TrackType.MUTATION)) {
            JMenuItem item = new JMenuItem("Sort by mutation count");
            item.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {

                    IGV.getInstance().sortByRegionScore(roi, RegionScoreType.MUTATION_COUNT, frame);
                    IGV.getInstance().getContentPane().repaint();

                }
            });

            add(item);
        }

        JMenuItem item = new JMenuItem("Sort by value");
        item.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {

                IGV.getInstance().sortByRegionScore(roi, RegionScoreType.SCORE, frame);
                IGV.getInstance().getContentPane().repaint();

            }
        });
        add(item);

        if (ScatterPlotUtils.hasPlottableTracks()) {
            addSeparator();
            item = new JMenuItem("Scatter Plot ...");
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {

                    String chr = roi.getChr();
                    int start = roi.getStart();
                    int end = roi.getEnd();
                    int zoom = frame.getZoom();
                    ScatterPlotUtils.openPlot(chr, start, end, zoom);
                }
            });
            add(item);
        }

    }


}
