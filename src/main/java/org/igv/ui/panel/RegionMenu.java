package org.igv.ui.panel;

import org.igv.charts.ScatterPlotUtils;
import org.igv.feature.RegionOfInterest;
import org.igv.track.RegionScoreType;
import org.igv.track.DataType;
import org.igv.ui.IGV;

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


        Set<DataType> loadedTypes = IGV.getInstance().getLoadedTypes();

        if (loadedTypes.contains(DataType.COPY_NUMBER) ||
                loadedTypes.contains(DataType.ALLELE_SPECIFIC_COPY_NUMBER) ||
                loadedTypes.contains(DataType.CNV)) {
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

        if (loadedTypes.contains(DataType.GENE_EXPRESSION)) {
            JMenuItem item = new JMenuItem("Sort by expression");
            item.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {

                    IGV.getInstance().sortByRegionScore(roi, RegionScoreType.EXPRESSION, frame);
                    IGV.getInstance().getContentPane().repaint();

                }
            });

            add(item);
        }

        if (loadedTypes.contains(DataType.MUTATION)) {
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
