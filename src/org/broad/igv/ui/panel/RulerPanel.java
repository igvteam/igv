/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/*
 * TrackPanel.java
 *
 * Created on Sep 5, 2007, 4:09:39 PM
 *
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.panel;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.NamedRunnable;

import javax.swing.*;
import javax.swing.event.MouseInputAdapter;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.geom.Rectangle2D;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 *         <p/>
 *         Lucida Blackletter
 *         Lucida Bright
 *         Lucida Calligraphy
 *         Lucida Fax
 *         Lucida Grande
 *         Lucida Handwriting
 *         Lucida Sans
 *         Lucida Sans Typewriter
 */
public class RulerPanel extends JPanel {

    private static Logger log = Logger.getLogger(RulerPanel.class);

    // TODO -- get from preferences
    boolean drawSpan = true;
    boolean drawEllipsis = false;
    private Font tickFont = FontManager.getFont(Font.BOLD, 9);
    private Font spanFont = FontManager.getFont(Font.BOLD, 12);
    private List<ClickLink> chromosomeRects = new ArrayList();

    private static Color dragColor = new Color(.5f, .5f, 1f, .3f);
    private static Color zoomBoundColor = new Color(0.5f, 0.5f, 0.5f);

    boolean dragging = false;
    int dragStart;
    int dragEnd;
    public static final String WHOLE_GENOME_TOOLTIP = "<html>Click on a chromosome number to jump to that chromosome," +
            "<br>or click and drag to zoom in.";
    public static final String CHROM_TOOLTIP = "Click and drag to zoom in.";

    ReferenceFrame frame;

    public RulerPanel(ReferenceFrame frame) {
        this.frame = frame;
        init();
    }

    private boolean isWholeGenomeView() {
        return getViewContext().getChrName().equals(Globals.CHR_ALL);
    }

    @Override
    protected void paintComponent(Graphics g) {

        super.paintComponent(g);

        g.setColor(Color.black);

        if (isWholeGenomeView()) {
            drawChromosomeTicks(g);
        } else {
            // Clear panel
            drawTicks(g);
            if (drawSpan) {
                drawSpan(g);
            }
            if (drawEllipsis) {
                drawEllipsis(g);
            }
        }


        if (dragging) {
            g.setColor(dragColor);
            int start = Math.min(dragStart, dragEnd);
            int w = Math.abs(dragEnd - dragStart);
            final int height = getHeight();
            g.fillRect(start, 0, w, height);

            g.setColor(zoomBoundColor);
            g.drawLine(dragStart, 0, dragStart, height);
            g.drawLine(dragEnd, 0, dragEnd, height);
        }

    }

    private ReferenceFrame getViewContext() {
        return frame;
    }

    private void drawSpan(Graphics g) {

        //TODO -- hack

        int w = getWidth();

        g.setFont(spanFont);

        // Positions are 1/2 open, subtract 1 to compensate
        int range = (int) (getViewContext().getScale() * w) + 1;

        TickSpacing ts = findSpacing(getViewContext(), range);
        String rangeString = formatNumber((double) range / ts.getUnitMultiplier()) + " " + ts.getMajorUnit();
        int strWidth = g.getFontMetrics().stringWidth(rangeString);
        int strHeight = g.getFontMetrics().getAscent();
        int strPosition = (w - strWidth) / 2;

        int lineY = getHeight() - 35 - strHeight / 2;

        g.drawLine(0, lineY, (w - strWidth) / 2 - 10, lineY);
        int[] arrowX = {0, 10, 10};
        int[] arrowY = {lineY, lineY + 3, lineY - 3};
        g.fillPolygon(arrowX, arrowY, arrowX.length);

        g.drawLine((w + strWidth) / 2 + 10, lineY, w, lineY);
        arrowX = new int[]{w, w - 10, w - 10};
        g.fillPolygon(arrowX, arrowY, arrowX.length);

        g.drawString(rangeString, strPosition, getHeight() - 35);

    }

    private void drawEllipsis(Graphics g) {
        double cytobandScale = ((double) getViewContext().getChromosomeLength()) / getWidth();

        double maxPixel = getViewContext().getMaxPixel();
        //visibleFraction = maxPixel < 0 ? 0 : ((double) getViewContext().getDataPanelWidth()) / maxPixel;

        int start = (int) ((getViewContext().getOrigin()) / cytobandScale);
        int span = (int) ((getWidth() * getViewContext().getScale()) / cytobandScale);
        int end = start + span;

        g.drawLine(start, 0, 0, getHeight());
        g.drawLine(end, 0, getWidth(), getHeight());

    }

    private void drawTicks(Graphics g) {

        int w = getWidth();
        if (w < 200) {
            return;
        }

        g.setFont(tickFont);

        int range = (int) (w * getViewContext().getScale());
        TickSpacing ts = findSpacing(getViewContext(), range);
        double spacing = ts.getMajorTick();

        // Find starting point closest to the current origin
        int nTick = (int) (getViewContext().getOrigin() / spacing) - 1;
        int l = (int) (nTick * spacing);
        int x = getViewContext().getScreenPosition(l - 1);    // 0 vs 1 based coordinates
        //int strEnd = Integer.MIN_VALUE;
        while (x < getWidth()) {
            l = (int) (nTick * spacing);
            x = getViewContext().getScreenPosition(l - 1);
            String chrPosition = formatNumber((double) l / ts.getUnitMultiplier()) +
                    " " + ts.getMajorUnit();
            int strWidth = g.getFontMetrics().stringWidth(chrPosition);
            int strPosition = x - strWidth / 2;
            //if (strPosition > strEnd) {

            if (nTick % 2 == 0) {
                g.drawString(chrPosition, strPosition, getHeight() - 15);
            }
            //strEnd = strPosition + strWidth;
            //}
            g.drawLine(x, getHeight() - 10, x, getHeight() - 2);
            nTick++;
        }
    }

    private void drawChromosomeTicks(Graphics g) {

        Font chrFont = FontManager.getFont(10);
        //this.removeAll();
        this.setLayout(null);

        // TODO -- remove hardcoded value
        int locationUnit = 1000;

        g.setFont(chrFont);

        Genome genome = IGV.getInstance().getGenomeManager().getCurrentGenome();
        if (genome == null) {
            log.info("No genome found with id: " + genome.getId());
            PreferenceManager.getInstance().remove(PreferenceManager.DEFAULT_GENOME_KEY);
        }

        boolean even = true;
        long offset = 0;
        chromosomeRects.clear();
        List<String> chrNames = genome.getChromosomeNames();
        if (chrNames == null) {
            log.info("No chromosomes found for genome: " + genome.getId());
            PreferenceManager.getInstance().remove(PreferenceManager.DEFAULT_GENOME_KEY);
        }
        if (chrNames.size() > 500) {
            return;
        }

        final FontMetrics fontMetrics = g.getFontMetrics();
        for (String chrName : chrNames) {
            Chromosome c = genome.getChromosome(chrName);
            int chrLength = c.getLength();

            int x = (int) (offset / (locationUnit * getViewContext().getScale()));
            int dw = (int) (chrLength / (locationUnit * getViewContext().getScale()));

            // Dont draw very small chromosome & contigs in whole genome view
            if (dw > 1) {

                g.drawLine(x, getHeight() - 10, x, getHeight() - 2);

                // Don't label chromosome if its width is < 5 pixels
                if (dw > 5) {
                    int center = x + dw / 2;

                    String displayName = null;
                    if (chrName.startsWith("gi|")) {
                        displayName = Genome.getNCBIName(chrName);
                    } else {
                        displayName = chrName.replace("chr", "");
                    }
                     int strWidth = fontMetrics.stringWidth(displayName);
                    int strPosition = center - strWidth / 2;


                    int y = (even ? getHeight() - 35 : getHeight() - 25);
                    g.drawString(displayName, strPosition, y);
                    int sw = (int) fontMetrics.getStringBounds(displayName, g).getWidth();
                    Rectangle clickRect = new Rectangle(strPosition, y - 15, sw, 15);
                    String tooltipText = "Jump to chromosome: " + chrName;
                    chromosomeRects.add(new ClickLink(clickRect, chrName, tooltipText));

                    even = !even;
                }
            }

            offset += chrLength;
        }
    }

    public static String formatNumber(double position) {

        //NumberFormatter f = new NumberFormatter();
        DecimalFormat formatter = new DecimalFormat();
        return formatter.format((int) position);
        //return f.valueToString(position);

    }

    public static TickSpacing findSpacing(ReferenceFrame frame, long maxValue) {

        if (maxValue < 10) {
            return new TickSpacing(1, "bp", 1);
        }
        // TODO -- hack, assumes location unit for whole genome is kilo-base
        boolean scaleInKB = frame.getChrName().equals(Globals.CHR_ALL);


        // Now man zeroes?
        int nZeroes = (int) Math.log10(maxValue);
        String majorUnit = scaleInKB ? "kb" : "bp";
        int unitMultiplier = 1;
        if (nZeroes > 9) {
            majorUnit = scaleInKB ? "tb" : "gb";
            unitMultiplier = 1000000000;
        }
        if (nZeroes > 6) {
            majorUnit = scaleInKB ? "gb" : "mb";
            unitMultiplier = 1000000;
        } else if (nZeroes > 3) {
            majorUnit = scaleInKB ? "mb" : "kb";
            unitMultiplier = 1000;
        }

        double nMajorTicks = maxValue / Math.pow(10, nZeroes - 1);
        if (nMajorTicks < 25) {
            return new TickSpacing(Math.pow(10, nZeroes - 1), majorUnit, unitMultiplier);
        } else {
            return new TickSpacing(Math.pow(10, nZeroes) / 2, majorUnit, unitMultiplier);
        }
    }

    private void init() {

        setBorder(BorderFactory.createLineBorder(UIConstants.TRACK_BORDER_GRAY));

        setCursor(Cursor.getDefaultCursor());
        if (isWholeGenomeView()) {
            this.setToolTipText(WHOLE_GENOME_TOOLTIP);
        } else {
            this.setToolTipText("Click and drag to zoom");
        }

        MouseInputAdapter mouseAdapter = new MouseInputAdapter() {

            int lastMousePressX;

            @Override
            public void mouseClicked(MouseEvent evt) {
                final MouseEvent e = evt;
                setCursor(Cursor.getDefaultCursor());                
                WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();
                try {

                    if (!isWholeGenomeView()) {
                        double newLocation = frame.getChromosomePosition(e.getX());
                        getViewContext().centerOnLocation(newLocation);
                    } else {
                        for (final ClickLink link : chromosomeRects) {
                            if (link.region.contains(e.getPoint())) {
                                NamedRunnable runnable = new NamedRunnable() {
                                    final String chrName = link.value;

                                    public void run() {

                                        getViewContext().setChrName(chrName);
                                        getViewContext().recordHistory();
                                        // TODO -- get rid of this ugly reference to IGV.theInstance
                                        IGV.getInstance().chromosomeChangeEvent(chrName);
                                    }

                                    public String getName() {
                                        return "Select chromosome: " + chrName;
                                    }
                                };

                                LongRunningTask.submit(runnable);

                                return;
                            }
                        }
                    }
                } finally {
                    WaitCursorManager.removeWaitCursor(token);
                }

            }

            @Override
            public void mouseMoved(MouseEvent e) {
                if (isWholeGenomeView()) {
                    for (ClickLink link : chromosomeRects) {
                        if (link.region.contains(e.getPoint())) {
                            setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
                            setToolTipText(link.tooltipText);
                            return;
                        }
                    }
                    setCursor(Cursor.getDefaultCursor());
                    setToolTipText(WHOLE_GENOME_TOOLTIP);

                }
            }

            @Override
            public void mouseEntered(MouseEvent e) {
                setCursor(Cursor.getDefaultCursor());
                if (isWholeGenomeView()) {
                    setToolTipText(WHOLE_GENOME_TOOLTIP);
                } else {
                    setToolTipText(CHROM_TOOLTIP);
                }

            }


            @Override
            public void mouseDragged(MouseEvent e) {
                if (Math.abs(e.getPoint().getX() - dragStart) > 1) {
                    dragEnd = e.getX();
                    dragging = true;
                    repaint();
                }
            }

            @Override
            public void mousePressed(MouseEvent e) {
                dragStart = e.getX();
            }

            @Override
            public void mouseReleased(MouseEvent e) {
                if (dragging) {
                    dragEnd = e.getX();
                    dragging = false;
                    zoom();
                }
            }

            //@Override
            //public void mouseExited(MouseEvent e) {
            //dragging = false;
            //}
        };

        addMouseMotionListener(mouseAdapter);

        addMouseListener(mouseAdapter);
    }

    private void zoom() {


        NamedRunnable runnable = new NamedRunnable() {
            public void run() {
                final ReferenceFrame vc = getViewContext();
                double s = vc.getChromosomePosition(dragStart);
                double e = vc.getChromosomePosition(dragEnd);
                if (e < s) {
                    double tmp = s;
                    s = e;
                    e = tmp;
                }
                if (e - s < 40) {
                    double c = (s + e) / 2;
                    s = c - 20;
                    e = c + 20;
                }

                String chr = null;
                if (isWholeGenomeView()) {
                    Genome genome = IGV.getInstance().getGenomeManager().getCurrentGenome();
                    Genome.ChromosomeCoordinate start = genome.getChromosomeCoordinate((int) s);
                    Genome.ChromosomeCoordinate end = genome.getChromosomeCoordinate((int) e);

                    chr = start.getChr();
                    s = start.getCoordinate();
                    e = end.getCoordinate();
                    if (end.getChr() != start.getChr()) {
                        e = genome.getChromosome(start.getChr()).getLength();
                    }
                } else {
                    chr = vc.getChrName();
                }

                vc.jumpTo(chr, (int) Math.min(s, e), (int) Math.max(s, e));
                vc.recordHistory();
            }

            public String getName() {
                return "Rule panel zoom";
            }
        };

        LongRunningTask.submit(runnable);

    }


    public static class TickSpacing {

        private double majorTick;
        private double minorTick;
        private String majorUnit = "";
        private int unitMultiplier = 1;

        TickSpacing(double majorTick, String majorUnit, int unitMultiplier) {
            this.majorTick = majorTick;
            this.minorTick = majorTick / 10.0;
            this.majorUnit = majorUnit;
            this.unitMultiplier = unitMultiplier;
        }

        public double getMajorTick() {
            return majorTick;
        }

        public double getMinorTick() {
            return minorTick;
        }

        public String getMajorUnit() {
            return majorUnit;
        }

        public void setMajorUnit(String majorUnit) {
            this.majorUnit = majorUnit;
        }

        public int getUnitMultiplier() {
            return unitMultiplier;
        }

        public void setUnitMultiplier(int unitMultiplier) {
            this.unitMultiplier = unitMultiplier;
        }
    }

// TODO -- possibly generalize?

    class ClickLink {

        Rectangle region;
        String value;
        String tooltipText;

        ClickLink(Rectangle region, String value, String tooltipText) {
            this.region = region;
            this.value = value;
            this.tooltipText = tooltipText;
        }
    }
}
