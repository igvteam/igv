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
 * TrackPanel.java
 *
 * Created on Sep 5, 2007, 4:09:39 PM
 *
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.panel;

import org.broad.igv.logging.*;
import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.ChromosomeCoordinate;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.sam.InsertionManager;
import org.broad.igv.sam.InsertionMarker;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.NamedRunnable;

import javax.swing.*;
import javax.swing.event.MouseInputAdapter;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.function.Consumer;

import static org.broad.igv.prefs.Constants.*;

/**
 * @author jrobinso
 */
public class RulerPanel extends JPanel {

    public static final Color INSERTION_MARKER_ZOOMEDOUT = new Color(64, 64, 64, 50);
    private static Logger log = LogManager.getLogger(RulerPanel.class);

    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat();

    private static final int INSERTION_ROW_HEIGHT = 9;

    // TODO -- get from preferences
    boolean drawSpan = true;
    private Font tickFont = FontManager.getFont(Font.BOLD, 9);
    private Font spanFont = FontManager.getFont(Font.BOLD, 12);

    private List<ClickLink> clickLinks = new ArrayList();

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
        return frame.getChrName().equals(Globals.CHR_ALL);
    }

    @Override
    protected void paintComponent(Graphics g) {

        super.paintComponent(g);

        if (PreferencesManager.getPreferences().getAntiAliasing()) {
            ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        }

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

        ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_DEFAULT);

    }

    private void render(Graphics g) {

        clickLinks.clear();
        g.setColor(Color.black);
        if (isWholeGenomeView()) {
            drawChromosomeTicks(g);
        } else {
            drawTicks(g);
            if (drawSpan) {
                drawSpan(g);
            }

            //  If zoomed in, use the 3rd gen visibility window
            Collection<AlignmentTrack> tracks = IGV.getInstance().getAlignmentTracks();
            int maxVizWindow = tracks.size() == 0 ? 0 : tracks.stream().mapToInt(t -> t.getVisibilityWindow()).max().getAsInt();
            if (!frame.getChrName().equals(Globals.CHR_ALL) && (frame.getEnd() - frame.getOrigin()) <= maxVizWindow) {
                drawInsertionMarkers(g);
            }
        }
    }

    private void drawSpan(Graphics g) {

        //TODO -- hack
        int w = getWidth();

        g.setFont(spanFont);

        int range = (int) (frame.getScale() * w) + 1;

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

    private void drawTicks(Graphics g) {

        int w = getWidth();
        if (w < 200) {
            return;
        }

        g.setFont(tickFont);

        // location unit for whole genome is kilobase
        boolean scaleInKB = frame.getChrName().equals(Globals.CHR_ALL);

        int range = (int) (w * frame.getScale());
        TickSpacing ts = findSpacing(range, scaleInKB);
        double spacing = ts.getMajorTick();

        // Find starting point closest to the current origin
        int nTick = (int) (frame.getExpansionAdjustedOrigin() / spacing) - 1;
        int l = (int) (nTick * spacing);
        int x = frame.getScreenPosition(l - 1 + 0.5);    // 0 vs 1 based coordinates, then center over base
        //int strEnd = Integer.MIN_VALUE;
        while (x < getWidth()) {
            l = (int) (nTick * spacing);
            x = frame.getScreenPosition(l - 1 + 0.5);
            String chrPosition = formatNumber((double) l / ts.getUnitMultiplier()) +
                    " " + ts.getMajorUnit();
            int strWidth = g.getFontMetrics().stringWidth(chrPosition);
            int strPosition = x - strWidth / 2;
            final int height = getHeight();

            if (nTick % 2 == 0 && strPosition > 0) {
                g.drawString(chrPosition, strPosition, height - 15);
            }
            if (x > 0) {
                g.drawLine(x, height - 10, x, height - 2);
            }

            nTick++;
        }
    }

    /**
     * Draw whole genome view, one "tick" per chromosome. *
     * @param g
     */
    private void drawChromosomeTicks(Graphics g) {

        Font chrFont = FontManager.getFont(10);
        this.setLayout(null);

        // Whole genome view uses kb units
        int locationUnit = 1000;

        g.setFont(chrFont);

        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        if (genome == null) {
            log.warn("No genome found");
            PreferencesManager.getPreferences().remove(DEFAULT_GENOME);
            return;
        }

        boolean even = true;
        List<String> chrNames = genome.getLongChromosomeNames();
        if (chrNames == null) {
            log.info("No chromosomes found for genome: " + PreferencesManager.getPreferences().getDefaultGenome());
            PreferencesManager.getPreferences().remove(DEFAULT_GENOME);
        }
        if (chrNames.size() > 500) {
            return;
        }

        final FontMetrics fontMetrics = g.getFontMetrics();
        for (final String chrName : chrNames) {
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

            if (x > 0) {
                g.drawLine(x, getHeight() - 10, x, getHeight() - 2);
            }

            // Don't label chromosome if its width is < 5 pixels
            if (dw > 5) {
                int center = x + dw / 2;

                String displayName = chrName;
                if (chrName.startsWith("gi|")) {
                    displayName = Genome.getNCBIName(chrName);
                } else if (chrName.length() < 6) {
                    displayName = chrName.replace("chr", "");
                }
                int strWidth = fontMetrics.stringWidth(displayName);
                int strPosition = center - strWidth / 2;


                int y = (even ? getHeight() - 35 : getHeight() - 25);
                g.drawString(displayName, strPosition, y);

                // Create clickable link if not dragging
                if (!dragging) {
                    int sw = (int) fontMetrics.getStringBounds(displayName, g).getWidth();
                    Rectangle clickRect = new Rectangle(strPosition, y - 15, sw, 15);
                    String tooltipText = "Jump to chromosome: " + chrName;
                    clickLinks.add(new ClickLink(clickRect, chrName, tooltipText, ch -> {
                        frame.changeChromosome((String) ch, true);
                    }));
                }

                even = !even;

            }
        }
    }

    private void drawInsertionMarkers(Graphics g) {

        List<InsertionMarker> insertionMarkers = InsertionManager.getInstance().getInsertions(frame.getChrName(), frame.getOrigin(), frame.getEnd() + 1);

        if (insertionMarkers == null) return;

        InsertionMarker expanded = frame.getExpandedInsertion();

        int w = (int) ((1.41 * INSERTION_ROW_HEIGHT) / 2);
        int y = getHeight() - INSERTION_ROW_HEIGHT - 1;

        for (InsertionMarker insertionMarker : insertionMarkers) {

            int x0;
            int x1;
            Polygon p;
            Color c;
            if (expanded != null && insertionMarker.position == expanded.position) {
                x0 = frame.getScreenPosition(insertionMarker.position);
                x1 = (int) ((insertionMarker.position + insertionMarker.size - frame.origin) / frame.getScale());
                p = new Polygon(
                        new int[]{x0 - w, x1 + w, x1, x0},
                        new int[]{y, y, y + INSERTION_ROW_HEIGHT, y + INSERTION_ROW_HEIGHT}, 4);
                c = Color.BLUE;

                String tooltipText = "Click to collapse insertion";
                clickLinks.add(new ClickLink(p, null, tooltipText, obj -> {
                    frame.setExpandedInsertion((InsertionMarker) obj);
                    IGV.getInstance().repaint();
                }));

            } else {
                x0 = frame.getScreenPosition(insertionMarker.position);
                x1 = x0;
                p = new Polygon(new int[]{x0 - w, x1 + w, x0},
                        new int[]{y, y, y + INSERTION_ROW_HEIGHT}, 3);

                double expandedInsertionWidth = insertionMarker.size  / frame.getScale();

                if(expandedInsertionWidth > 5) {
                    c = Color.BLACK ;
                    String tooltipText = "Click to expand insertion (" + (insertionMarker.size + "bases)");
                    Rectangle clickArea = p.getBounds();
                    clickArea.y -= 2;
                    clickArea.height += 2;
                    clickLinks.add(new ClickLink(clickArea, insertionMarker, tooltipText, obj -> {
                        frame.setExpandedInsertion((InsertionMarker) obj);
                        IGV.getInstance().repaint();
                    }));
                } else {
                    c = INSERTION_MARKER_ZOOMEDOUT;
                }
            }

            if (x1 + w >= 0 && x0 - w <= getWidth()) {
                g.setColor(c);
                g.fillPolygon(p);
            }
        }
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

        //setBorder(BorderFactory.createLineBorder(UIConstants.TRACK_BORDER_GRAY));
        setCursor(Cursor.getDefaultCursor());
        if (isWholeGenomeView()) {
            this.setToolTipText(WHOLE_GENOME_TOOLTIP);
        } else {
            this.setToolTipText("Click and drag to zoom");
        }

        MouseInputAdapter mouseAdapter = new MouseInputAdapter() {

            private MouseEvent mouseDown;

            @Override
            public void mouseMoved(MouseEvent e) {
                for (ClickLink link : clickLinks) {
                    if (link.region.contains(e.getPoint())) {
                        setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
                        setToolTipText(link.tooltipText);
                        return;
                    }
                }
                setCursor(Cursor.getDefaultCursor());
                setToolTipText(isWholeGenomeView() ? WHOLE_GENOME_TOOLTIP : CHROM_TOOLTIP);

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
            public void mousePressed(final MouseEvent e) {
                if (e.isPopupTrigger()) {
                    // ignore
                } else {
                    dragStart = e.getX();
                    mouseDown = e;
                }
            }


            @Override
            public void mouseReleased(MouseEvent e) {
                if (e.isPopupTrigger()) {
                    // ignore
                } else {
                    if (dragging) {
                        dragEnd = e.getX();
                        dragging = false;
                        zoom();
                    } else {
                        if (mouseDown != null && distance(mouseDown, e) < 5) {
                            doMouseClick(e);
                        }
                    }
                }
            }

            private double distance(MouseEvent e1, MouseEvent e2) {
                double dx = e1.getX() - e2.getX();
                double dy = e1.getY() - e2.getY();
                return Math.sqrt(dx * dx + dy * dy);
            }


            /**
             * mouseClick is not used because in Java a mouseClick is emitted ONLY if the mouse has not moved
             * between press and release.  This is really difficult to achieve, even if trying.
             *
             * @param evt
             */
            public void doMouseClick(MouseEvent evt) {
                final MouseEvent e = evt;
                setCursor(Cursor.getDefaultCursor());
                WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();
                try {
                    boolean clickHandled = false;
                    for (final ClickLink link : clickLinks) {
                        if (link.region.contains(e.getPoint())) {
                            link.action.accept(link.value);
                            clickHandled = true;
                        }
                    }
                    if (!clickHandled && !isWholeGenomeView()) {
                        double newLocation = frame.getChromosomePosition(e);
                        frame.centerOnLocation(newLocation);
                    }

                } finally {
                    WaitCursorManager.removeWaitCursor(token);
                }

            }

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
                Genome genome = GenomeManager.getInstance().getCurrentGenome();

                if (isWholeGenomeView()) {

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
                    s = Math.max(0, s);
                    e = Math.min(genome.getChromosome(chr).getLength(), e);

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

        Shape region;
        Object value;
        String tooltipText;

        Consumer action;

        ClickLink(Shape region, Object value, String tooltipText, Consumer action) {
            this.region = region;
            this.value = value;
            this.tooltipText = tooltipText;
            this.action = action;
        }

    }


}
