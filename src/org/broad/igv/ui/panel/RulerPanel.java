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
import org.broad.igv.dev.affective.AffectiveUtils;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.exome.ExomeBlock;
import org.broad.igv.feature.exome.ExomeReferenceFrame;
import org.broad.igv.feature.genome.ChromosomeCoordinate;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeImpl;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.renderer.GraphicUtils;
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
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * @author jrobinso
 */
public class RulerPanel extends JPanel {

    private static Logger log = Logger.getLogger(RulerPanel.class);

    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat();

    private static Color grey1 = new Color(120, 120, 120);
    private static Color grey2 = new Color(200, 200, 200);
    private static Color gene1 = new Color(0, 0, 150, 150);
    private static Color gene2 = new Color(0, 150, 0, 150);

    // TODO -- get from preferences
    boolean affective = false;
    boolean drawSpan = true;
    boolean drawEllipsis = false;
    private Font tickFont = FontManager.getFont(Font.BOLD, 9);
    private Font spanFont = FontManager.getFont(Font.BOLD, 12);

    private List<ClickLink> chromosomeRects = new ArrayList();
    private List<MouseRect> mouseRects = new ArrayList<MouseRect>();


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
        affective = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.AFFECTIVE_ENABLE);
        drawSpan = !affective;
        init();
    }

    private boolean isWholeGenomeView() {
        return frame.getChrName().equals(Globals.CHR_ALL);
    }

    @Override
    protected void paintComponent(Graphics g) {

        super.paintComponent(g);

        ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        render(g);


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

    private void render(Graphics g) {
        g.setColor(Color.black);

        if (isWholeGenomeView()) {
            drawChromosomeTicks(g);
        } else if (FrameManager.isExomeMode()) {
            // TODO -- hack,
            ExomeReferenceFrame exomeFrame = (ExomeReferenceFrame) FrameManager.getDefaultFrame();
            drawExomeBlocks(g, exomeFrame);
            drawExomeGenes(g, exomeFrame);
            if (drawSpan) {
                drawSpan(g);
            }

        } else {
            // Clear panel
            if (affective) {
                drawTimeTicks(g);
            } else {
                drawTicks(g);
            }
            if (drawSpan) {
                drawSpan(g);
            }
            if (drawEllipsis) {
                drawEllipsis(g);
            }
        }

    }

    private void drawSpan(Graphics g) {

        //TODO -- hack
        int w = getWidth();

        g.setFont(spanFont);


        int range;
        if (FrameManager.isExomeMode()) {
            range = (int) (frame.getEnd() - frame.getOrigin()) + 1;
        } else {
            range = (int) (frame.getScale() * w) + 1;
        }

        // TODO -- hack, assumes location unit for whole genome is kilo-base
        boolean scaleInKB = frame.getChrName().equals(Globals.CHR_ALL);

        TickSpacing ts = findSpacing(range, scaleInKB);
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
        double cytobandScale = ((double) frame.getChromosomeLength()) / getWidth();

        double maxPixel = frame.getMaxPixel();
        //visibleFraction = maxPixel < 0 ? 0 : ((double) getViewContext().getDataPanelWidth()) / maxPixel;

        int start = (int) ((frame.getOrigin()) / cytobandScale);
        int span = (int) ((getWidth() * frame.getScale()) / cytobandScale);
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

        // TODO -- hack, assumes location unit for whole genome is kilo-base
        boolean scaleInKB = frame.getChrName().equals(Globals.CHR_ALL);

        int range = (int) (w * frame.getScale());
        TickSpacing ts = findSpacing(range, scaleInKB);
        double spacing = ts.getMajorTick();

        // Find starting point closest to the current origin
        int nTick = (int) (frame.getOrigin() / spacing) - 1;
        int l = (int) (nTick * spacing);
        int x = frame.getScreenPosition(l - 1);    // 0 vs 1 based coordinates
        //int strEnd = Integer.MIN_VALUE;
        while (x < getWidth()) {
            l = (int) (nTick * spacing);
            x = frame.getScreenPosition(l - 1);
            String chrPosition = formatNumber((double) l / ts.getUnitMultiplier()) +
                    " " + ts.getMajorUnit();
            int strWidth = g.getFontMetrics().stringWidth(chrPosition);
            int strPosition = x - strWidth / 2;
            //if (strPosition > strEnd) {

            final int height = getHeight();
            if (nTick % 2 == 0) {
                g.drawString(chrPosition, strPosition, height - 15);
            }

            g.drawLine(x, height - 10, x, height - 2);
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

        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        if (genome == null) {
            log.info("No genome found with id: " + genome.getId());
            PreferenceManager.getInstance().remove(PreferenceManager.DEFAULT_GENOME_KEY);
        }

        boolean even = true;
        long offset = 0;
        chromosomeRects.clear();
        List<String> chrNames = genome.getLongChromosomeNames();
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
            if (c == null) {
                log.info("Chromosome '" + chrName + "' not found");
                continue;
            }
            int chrLength = c.getLength();

            double scale = frame.getScale();
            int gStart = genome.getGenomeCoordinate(chrName, 0);
            int x = (int) (gStart / scale);
            int dw = (int) (chrLength / (locationUnit * scale));


            g.drawLine(x, getHeight() - 10, x, getHeight() - 2);

            // Don't label chromosome if its width is < 5 pixels
            if (dw > 5) {
                int center = x + dw / 2;

                String displayName = null;
                if (chrName.startsWith("gi|")) {
                    displayName = GenomeImpl.getNCBIName(chrName);
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

            offset += chrLength;
        }
    }

    private void drawExomeBlocks(Graphics g, ExomeReferenceFrame frame) {

        String chr = frame.getChrName();
        List<ExomeBlock> blocks = frame.getBlocks(chr);
        int idx = frame.getFirstBlockIdx();
        ExomeBlock b;

        Rectangle visibleRect = this.getVisibleRect();


        int lastPStart = -1;
        int pStart;
        int pEnd;
        int exomeOrigin = frame.getExomeOrigin();
        int visibleBlockCount = 0;
        int blockGap = 0; // TODO --
        int top = visibleRect.y + visibleRect.height - 10;
        do {
            b = blocks.get(idx);

            pStart = (int) ((b.getExomeStart() - exomeOrigin) / frame.getScale()) + visibleBlockCount * blockGap;
            pEnd = (int) ((b.getExomeEnd() - exomeOrigin) / frame.getScale()) + visibleBlockCount * blockGap;

            // Don't draw over previously drawn region -- can happen when zoomed out.
            if (pEnd > lastPStart) {

                lastPStart = pStart;
                if (pEnd == pStart) pEnd++;


                b.setScreenBounds(pStart, pEnd);

                Rectangle rect = new Rectangle(pStart, top, pEnd - pStart, 10);


                Graphics2D exomeGraphics = (Graphics2D) g.create();
                exomeGraphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

//Shape clip = exomeGraphics.getClip();

                Color c = idx % 2 == 0 ? grey1 : grey2;

                exomeGraphics.setColor(c);
                exomeGraphics.fill(rect);
                //GraphicUtils.drawCenteredText(String.valueOf(idx), rect, exomeGraphics);


                visibleBlockCount++;
            }
            idx++;


        }
        while ((pStart < visibleRect.x + visibleRect.width) && idx < blocks.size());


    }

    private void drawExomeGenes(Graphics g, ExomeReferenceFrame frame) {

        mouseRects.clear();

        String chr = frame.getChrName();
        List<ExomeReferenceFrame.Gene> genes = frame.getGenes(chr);

        int idx = FeatureUtils.getIndexBefore(frame.getOrigin(), genes);
        Rectangle visibleRect = this.getVisibleRect();
        FontMetrics fm = g.getFontMetrics();

        int lastPStart = -1;
        int pStart;
        int pEnd;
        int exomeOrigin = frame.getExomeOrigin();
        int visibleBlockCount = 0;
        int blockGap = 0; // TODO --
        int top = visibleRect.y + visibleRect.height - 30;
        do {
            ExomeReferenceFrame.Gene gene = genes.get(idx);

            double exomeStart = frame.genomeToExomePosition(gene.getStart());
            double exomeEnd = frame.genomeToExomePosition(gene.getEnd());

            pStart = (int) ((exomeStart - exomeOrigin) / frame.getScale()) + visibleBlockCount * blockGap;
            pEnd = (int) ((exomeEnd - exomeOrigin) / frame.getScale()) + visibleBlockCount * blockGap;

            // Don't draw over previously drawn region -- can happen when zoomed out.
            if (pEnd > lastPStart) {

                lastPStart = pStart;
                if (pEnd == pStart) pEnd++;

                Rectangle rect = new Rectangle(pStart, top, pEnd - pStart, 20);


                Graphics2D exomeGraphics = (Graphics2D) g.create();
                exomeGraphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
//Shape clip = exomeGraphics.getClip();

                Color c = idx % 2 == 0 ? gene1 : gene2;

                exomeGraphics.setColor(c);
                exomeGraphics.fill(rect);

                Rectangle2D sb = fm.getStringBounds(gene.getName(), g);
                if (sb.getWidth() < rect.getWidth()) {
                    exomeGraphics.setColor(Color.black);
                    GraphicUtils.drawCenteredText(gene.getName(), rect, exomeGraphics);
                }

                mouseRects.add(new MouseRect(rect, gene.getName()));

                visibleBlockCount++;
            }
            idx++;


        }
        while ((pStart < visibleRect.x + visibleRect.width) && idx < genes.size());

        Collections.sort(mouseRects, new Comparator<MouseRect>() {
            @Override
            public int compare(MouseRect mr1, MouseRect mr2) {
                return mr1.width() - mr2.width();
            }
        });


    }

    public static String formatNumber(double position) {

        //NumberFormatter f = new NumberFormatter();
        return DECIMAL_FORMAT.format((int) position);
        //return f.valueToString(position);

    }

    public static TickSpacing findSpacing(long maxValue, boolean scaleInKB) {

        if (maxValue < 10) {
            return new TickSpacing(1, "bp", 1);
        }


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
                        frame.centerOnLocation(newLocation);
                    } else {
                        for (final ClickLink link : chromosomeRects) {
                            if (link.region.contains(e.getPoint())) {
                                NamedRunnable runnable = new NamedRunnable() {
                                    final String chrName = link.value;

                                    public void run() {

                                        frame.setChromosomeName(chrName);
                                        frame.recordHistory();
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

                } else {
                    for (MouseRect mr : mouseRects) {
                        if (mr.containsPoint(e.getPoint())) {
                            RulerPanel.this.setToolTipText(mr.getText());
                            return;
                        }
                    }
                    RulerPanel.this.setToolTipText(CHROM_TOOLTIP);

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
                double s = frame.getChromosomePosition(dragStart);
                double e = frame.getChromosomePosition(dragEnd);
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

                s = Math.max(0.0, s);
                String chr = null;
                if (isWholeGenomeView()) {
                    Genome genome = GenomeManager.getInstance().getCurrentGenome();
                    ChromosomeCoordinate start = genome.getChromosomeCoordinate((int) s);
                    ChromosomeCoordinate end = genome.getChromosomeCoordinate((int) e);

                    chr = start.getChr();
                    s = start.getCoordinate();
                    e = end.getCoordinate();
                    if (end.getChr() != start.getChr()) {
                        e = genome.getChromosome(start.getChr()).getLength();
                    }
                } else {
                    chr = frame.getChrName();
                }

                frame.jumpTo(chr, (int) Math.min(s, e), (int) Math.max(s, e));
                frame.recordHistory();
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


    /**
     * Special renderer for "Affective Computing" timescale,  chromosome => day, units are hours, minutes, seconds.
     *
     * @param g
     */
    private void drawTimeTicks(Graphics g) {

        double timeStep = 1.0 / AffectiveUtils.POINTS_PER_SECOND;

        int w = getWidth();
        double start = frame.getOrigin();
        double end = frame.getEnd();
        double seconds = (start - end) * timeStep;

        // Determine step sizes
        double secsPerPixel = frame.getScale() * timeStep;
        double minsPerPixel = secsPerPixel / 60;
        double hoursPerPixel = minsPerPixel / 60;

        g.setFont(tickFont);
        FontMetrics fm = g.getFontMetrics();

        int startHour = (int) ((start * timeStep) / 3600);
        int endHour = (int) ((end * timeStep) / 3600) + 1;

        double originHour = (start * timeStep) / 3600;


        // Rectangle rect = getBounds();

        int height = getHeight();

        for (double h = startHour; h < endHour; h++) {
            double pixel = (int) ((h - originHour) / hoursPerPixel);

            if (pixel > w) {
                break;
            }
            if (pixel > 0) {
                g.drawLine((int) pixel, height, (int) pixel, height - 15);

                // Label
                int absoluteHour = AffectiveUtils.START_TIME_HR + (int) h;
                String label = absoluteHour + ":00";
                int labelWidth = fm.stringWidth(label);
                int labelX = (int) pixel - labelWidth / 2;
                if (labelX > 0) {
                    g.drawString(label, labelX, height - 20);
                }
            }

            // If room for 1/4 hours
            if (15 / minsPerPixel > 2) {
                pixel = (int) ((h - originHour) / hoursPerPixel);
                for (int mm = 0; mm < 60; mm += 15) {
                    double dx = mm / minsPerPixel;
                    int mPixel = (int) (pixel + dx);
                    if (mPixel > w) {
                        break;
                    }
                    if (mPixel > 0) {
                        g.drawLine(mPixel, height, mPixel, height - 10);
                    }

                }
            }

            // If room for minutes do minutes
            pixel = (int) ((h - originHour) / hoursPerPixel);
            if (1 / minsPerPixel > 4) {
                for (int m = 1; m < 60; m++) {
                    double dx = m / minsPerPixel;
                    int mPixel = (int) (pixel + dx);
                    if (mPixel > w) {
                        break;
                    }
                    if (mPixel > 0) {
                        g.drawLine(mPixel, height, mPixel, height - 5);
                    }
                }

                // Seconds
                if (1 / secsPerPixel > 4) {
                    for (int s = 1; s < 60; s++) {
                        double dx = s / secsPerPixel;
                        int sPixel = (int) (pixel + dx);
                        if (sPixel > w) {
                            break;
                        }
                        if (sPixel > 0) {
                            g.drawLine(sPixel, height, sPixel, height - 5);
                        }
                    }
                }
            }


        }
    }

    static class MouseRect {
        Rectangle bounds;
        String text;

        MouseRect(Rectangle bounds, String text) {
            this.bounds = bounds;
            this.text = text;
        }

        boolean containsPoint(Point p) {
            return bounds.contains(p);
        }

        String getText() {
            return text;
        }

        public int width() {
            return bounds.width;
        }
    }


}
