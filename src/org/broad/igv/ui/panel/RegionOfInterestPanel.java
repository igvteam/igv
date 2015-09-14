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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.panel;

import org.broad.igv.Globals;
import org.broad.igv.util.blat.BlatClient;
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
        popupMenu.add(item);
        // Disable copySequence if region exceeds 1 MB
        if (roi.getEnd() - roi.getStart() > 1000000) {
            item.setEnabled(false);
        }
        popupMenu.add(item);

        item = new JMenuItem("Blat sequence");
        item.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                BlatClient.doBlatQuery(roi.getChr(), roi.getStart(), roi.getEnd());
            }
        });
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
