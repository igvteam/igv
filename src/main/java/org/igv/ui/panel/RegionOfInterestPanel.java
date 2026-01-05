package org.igv.ui.panel;

import org.igv.Globals;
import org.igv.feature.Strand;
import org.igv.sam.AlignmentInterval;
import org.igv.sam.AlignmentRenderer;
import org.igv.sam.InsertionManager;
import org.igv.sam.InsertionMarker;
import org.igv.track.RenderContext;
import org.igv.util.blat.BlatClient;
import org.igv.feature.RegionOfInterest;
import org.igv.feature.genome.Genome;
import org.igv.feature.genome.GenomeManager;
import org.igv.ui.IGV;
import org.igv.util.LongRunningTask;
import org.igv.util.NamedRunnable;

import javax.swing.*;
import javax.swing.event.MouseInputAdapter;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.util.Collection;
import java.util.List;

/**
 * @author eflakes
 */
public class RegionOfInterestPanel extends JPanel {

    private final boolean darkMode;
    ReferenceFrame frame;
    RegionOfInterest focusROI;
    boolean switchStartOrEnd;

    // There can only be 1 selected region, irrespective of the number of panels
    private static RegionOfInterest selectedRegion = null;

    public RegionOfInterestPanel(ReferenceFrame frame) {
        this.darkMode = Globals.isDarkMode();
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

        //drawInsertionMarkers((Graphics2D) g, getHeight());

        g.setColor(darkMode ? Color.WHITE : Color.BLACK);
        g.drawLine(0, 0, getWidth(), 0);
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
                        IGV.copySequenceToClipboard(genome, roi.getChr(), roi.getStart(), roi.getEnd(), Strand.NONE);
                    }
                });
            }
        });
        popupMenu.add(item);
        // Disable copySequence if region exceeds 1 MB
        final int roiLength = roi.getEnd() - roi.getStart();
        if (roiLength > 1000000) {
            item.setEnabled(false);
        }
        popupMenu.add(item);

        item = new JMenuItem("BLAT sequence");
        if (roiLength > 20 && roiLength < 8000) {
            item.addActionListener(e -> BlatClient.doBlatQueryFromRegion(roi.getChr(), roi.getStart(), roi.getEnd(), Strand.NONE));
        } else {
            item.setEnabled(false);
        }
        popupMenu.add(item);


        popupMenu.add(new JSeparator());

        item = new JMenuItem("Delete");
        item.addActionListener(e -> {
            IGV.getInstance().getSession().getRegionsOfInterest(frame.getChrName()).remove(roi);
            IGV.getInstance().repaint();
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

        boolean dragging = false;

        @Override
        public void mousePressed(MouseEvent e) {
            if ((e.getModifiers() & MouseEvent.CTRL_MASK) != 0) {
                focusROI = getRegionOfInterest(e.getX());
                if (focusROI != null) {
                    int curPos = (int) frame.getChromosomePosition(e);
                    int startDist = Math.abs(focusROI.getStart() - curPos);
                    int endDist = Math.abs(focusROI.getEnd() - curPos);
                    if (startDist < endDist) {
                        switchStartOrEnd = true;
                    } else {
                        switchStartOrEnd = false;
                    }
                }
            } else {
                showPopup(e);
            }
        }

        @Override
        public void mouseDragged(MouseEvent e) {
            dragging = true;
            if (focusROI != null) {
                if (switchStartOrEnd) {
                    focusROI.setStart((int) frame.getChromosomePosition(e));
                } else {
                    focusROI.setEnd((int) frame.getChromosomePosition(e));
                }
                IGV.getInstance().repaint();
            }
        }

        @Override
        public void mouseReleased(MouseEvent e) {
            if (dragging && selectedRegion != null) {
                selectedRegion = null;
                IGV.getInstance().repaint();
            }
            focusROI = null;
            dragging = false;
        }

        @Override
        public void mouseMoved(MouseEvent e) {
            RegionOfInterest roi = getRegionOfInterest(e.getX());
            if (roi != null) {
                setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
                setToolTipText("<html>" + roi.getTooltip() + "<br>To resize use ctrl-click-drag.");
                if (selectedRegion != roi) {
                    selectedRegion = roi;
                    IGV.getInstance().repaint();
                }

            } else {
                if (selectedRegion != null) {
                    selectedRegion = null;
                    IGV.getInstance().repaint();
                }
                setToolTipText("");
                setCursor(Cursor.getDefaultCursor());
            }
        }

        @Override
        public void mouseExited(MouseEvent mouseEvent) {
            if (!dragging && selectedRegion != null) {
                selectedRegion = null;
                IGV.getInstance().repaint();
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
