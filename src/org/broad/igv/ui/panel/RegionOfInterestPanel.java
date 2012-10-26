/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.panel;

import org.broad.igv.Globals;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.NamedRunnable;

import javax.swing.*;
import javax.swing.event.MouseInputAdapter;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.util.Collection;

/**
 * @author eflakes
 */
public class RegionOfInterestPanel extends JPanel {

    PopupMenu popup;

    ReferenceFrame frame;

    // There can only be 1 selected region, irrespective of the number of panels
    private static RegionOfInterest selectedRegion = null;

    public RegionOfInterestPanel(ReferenceFrame frame) {

        setToolTipText("Regions of Interest");
        this.frame = frame;
        MouseInputAdapter ma = new ROIMouseAdapater();
        addMouseListener(ma);
        addMouseMotionListener(ma);
    }


    @Override
    public void paintComponent(final Graphics g) {

        super.paintComponent(g);
        ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        // Draw regions of interest?
        drawRegionsOfInterest((Graphics2D) g, getHeight());

        g.setColor(Color.BLACK);
        g.drawRect(0, 0, getWidth(), getHeight());
    }


    public void drawRegionsOfInterest(final Graphics2D g, int height) {

        Collection<RegionOfInterest> regions = getRegions();

        if (regions == null || regions.isEmpty()) {
            return;
        }


        for (RegionOfInterest regionOfInterest : regions) {

            int regionStart = regionOfInterest.getStart();
            int regionEnd = regionOfInterest.getEnd();

            // This is ugly, but neccessary the way the "whole genome" is treated as another chromosome
            if (frame.getChrName().equals(Globals.CHR_ALL)) {
                Genome genome = GenomeManager.getInstance().getCurrentGenome();
                regionStart = genome.getGenomeCoordinate(regionOfInterest.getChr(), regionStart);
                regionEnd = genome.getGenomeCoordinate(regionOfInterest.getChr(), regionEnd);
            }

            int start = frame.getScreenPosition(regionStart);
            int end = frame.getScreenPosition(regionEnd);
            int regionWidth = Math.max(1, end - start);

            g.setColor(regionOfInterest.getBackgroundColor());
            g.fillRect(start, 0, regionWidth, height);

        }
    }


    /**
     * Return the region of interest at the screen pixel location.
     *
     * @param px
     * @return
     */
    RegionOfInterest getRegionOfInterest(int px) {

        double pos = frame.getChromosomePosition(px);

        Collection<RegionOfInterest> roiList = getRegions();
        if (roiList != null) {
            for (RegionOfInterest roi : roiList) {
                if (pos > roi.getStart() && pos < roi.getEnd()) {
                    return roi;
                }
            }
        }

        return null;

    }

    protected static JPopupMenu getPopupMenu(final Component parent, final RegionOfInterest roi, final ReferenceFrame frame) {

        //Set<TrackType> loadedTypes = IGV.getInstance().getLoadedTypes();

        JPopupMenu popupMenu = new RegionMenu(roi, frame);

        popupMenu.addSeparator();

        JMenuItem item = new JMenuItem("Zoom");
        item.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {

                frame.jumpTo(roi.getChr(), roi.getStart(), roi.getEnd());

                String locusString = roi.getLocusString();
                IGV.getInstance().getSession().getHistory().push(locusString, frame.getZoom());

            }
        });
        popupMenu.add(item);

        item = new JMenuItem("Edit description...");
        item.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                String desc = JOptionPane.showInputDialog(parent, "Add or edit region description:", roi.getDescription());
                roi.setDescription(desc);
                IGV.getInstance().getSession().getRegionsOfInterestObservable().setChangedAndNotify();

            }
        });
        popupMenu.add(item);

        item = new JMenuItem("Copy sequence");
        item.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {

                LongRunningTask.submit(new NamedRunnable() {
                    public String getName() {
                        return "Copy sequence";
                    }

                    public void run() {
                        Genome genome = GenomeManager.getInstance().getCurrentGenome();
                        IGV.copySequenceToClipboard(genome, roi.getChr(), roi.getStart(), roi.getEnd());
                    }
                });


            }
        });

        // Disable copySequence if region exceeds a MB
        if (roi.getEnd() - roi.getStart() > 1000000) {
            item.setEnabled(false);
        }
        popupMenu.add(item);


        popupMenu.add(new JSeparator());

        item = new JMenuItem("Delete");
        item.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                IGV.getInstance().getSession().getRegionsOfInterest(frame.getChrName()).remove(roi);
                IGV.getInstance().repaintDataAndHeaderPanels();
            }
        });
        popupMenu.add(item);

        return popupMenu;
    }


    public static RegionOfInterest getSelectedRegion() {
        return selectedRegion;
    }

    public static void setSelectedRegion(RegionOfInterest region) {
        selectedRegion = region;
    }

    class ROIMouseAdapater extends MouseInputAdapter {

        @Override
        public void mousePressed(MouseEvent e) {
            showPopup(e);
        }


        @Override
        public void mouseReleased(MouseEvent e) {
            //showPopup(e);
        }


        @Override
        public void mouseMoved(MouseEvent e) {
            RegionOfInterest roi = getRegionOfInterest(e.getX());
            if (roi != null) {
                setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
                setToolTipText(roi.getTooltip());
                if (selectedRegion != roi) {
                    selectedRegion = roi;
                    IGV.getInstance().repaintDataPanels();
                }

            } else {
                if (selectedRegion != null) {
                    selectedRegion = null;
                    IGV.getInstance().repaintDataPanels();
                }
                setToolTipText("");
                setCursor(Cursor.getDefaultCursor());
            }
        }

        @Override
        public void mouseExited(MouseEvent mouseEvent) {
            if (selectedRegion != null) {
                selectedRegion = null;
                IGV.getInstance().repaintDataPanels();
            }
        }

        private void showPopup(MouseEvent e) {

            RegionOfInterest roi = getRegionOfInterest(e.getX());
            if (roi != null) {

                getPopupMenu(RegionOfInterestPanel.this, roi, frame).show(e.getComponent(), e.getX(), e.getY());
            }

        }
    }


    /**
     * A convenience method for returning the regions of interest for the current frame.
     */

    private Collection<RegionOfInterest> getRegions() {
        return IGV.getInstance().getSession().getRegionsOfInterest(frame.getChrName());
    }
}
