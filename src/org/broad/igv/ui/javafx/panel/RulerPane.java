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

package org.broad.igv.ui.javafx.panel;

import javafx.geometry.Bounds;
import javafx.geometry.Point2D;
import javafx.geometry.Rectangle2D;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.text.Text;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sam.InsertionManager;
import org.broad.igv.sam.InsertionMarker;
import org.broad.igv.ui.javafx.FontMetrics;
import org.broad.igv.ui.javafx.IGVBackendPlaceholder;
import org.broad.igv.ui.javafx.ResizableCanvas;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import static org.broad.igv.prefs.Constants.DEFAULT_GENOME;

/**
 * @author jrobinso
 * @author eby JavaFX port
 */
// TODO: Not dealing with DnD for now.
public class RulerPane extends ResizableCanvas {

    private static Logger log = Logger.getLogger(RulerPane.class);

    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat();

    // TODO -- get from preferences
    boolean drawSpan = true;
    private Font tickFont = IGVBackendPlaceholder.getFont(FontWeight.BOLD, 9);
    private Font spanFont = IGVBackendPlaceholder.getFont(FontWeight.BOLD, 12);
    private Font chrFont = IGVBackendPlaceholder.getFont(10);

    private List<ClickLink> chromosomeRects = new ArrayList<ClickLink>();
    private List<MouseRect> mouseRects = new ArrayList<MouseRect>();


    private static Color dragColor = Color.color(.5f, .5f, 1f, .3f);
    private static Color zoomBoundColor = Color.color(0.5f, 0.5f, 0.5f);

    boolean dragging = false;
    int dragStart;
    int dragEnd;
    public static final String WHOLE_GENOME_TOOLTIP = "<html>Click on a chromosome number to jump to that chromosome," +
            "<br>or click and drag to zoom in.";
    public static final String CHROM_TOOLTIP = "Click and drag to zoom in.";

    ReferenceFrame frame;

    public RulerPane(ReferenceFrame frame) {
        this.frame = frame;
        setMinHeight(80);
        setMaxHeight(80);
        setPrefHeight(80);
        init();

        // Re-render on change of width/height or chromosome.  Note that the height is fixed and
        // so that listener should never execute.  However, leaving this in place as a pattern
        // and also because it may not be fixed in the long run.
        frame.chromosomeNameProperty().addListener((observable, oldValue, newValue) -> render());
        frame.zoomProperty().addListener((observable, oldValue, newValue) -> render());
        this.prefWidthProperty().addListener((observable, oldValue, newValue) -> render());
        this.prefHeightProperty().addListener((observable, oldValue, newValue) -> render());

        render();
    }

    private boolean isWholeGenomeView() {
        return frame.getChrName().equals(Globals.CHR_ALL);
    }

    public void render() {
        log.info("rendering chr: " + frame.getChrName());

        Canvas canvas = getCanvas();
        GraphicsContext graphicsContext = canvas.getGraphicsContext2D();
        graphicsContext.clearRect(0.0, 0.0, canvas.getWidth(), canvas.getHeight());

        graphicsContext.setStroke(Color.rgb(200, 200, 210));
        graphicsContext.strokeRect(0.0, 0.0, canvas.getWidth(), canvas.getHeight());

        graphicsContext.setFill(Color.BLACK);
        graphicsContext.setStroke(Color.BLACK);

        if (isWholeGenomeView()) {
            drawChromosomeTicks(graphicsContext);
        } else {

            InsertionMarker i = InsertionManager.getInstance().getSelectedInsertion(frame.getChrName());

            drawTicks(graphicsContext, i);

            if (drawSpan) {
                drawSpan(graphicsContext, i);
            }
        }
    }

    private void drawSpan(GraphicsContext graphicsContext, InsertionMarker i) {

        //TODO -- hack
        double w = getPrefWidth();

        graphicsContext.setFont(spanFont);

        double range = (frame.getScale() * w) + 1;

        // TODO -- hack, assumes location unit for whole genome is kilo-base
        boolean scaleInKB = frame.getChrName().equals(Globals.CHR_ALL);

        TickSpacing ts = findSpacing(range, scaleInKB);
        String rangeString = formatNumber((double) range / ts.getUnitMultiplier()) + " " + ts.getMajorUnit();
        Bounds bounds = FontMetrics.getBoundsInFont(rangeString, spanFont);
        double strWidth = bounds.getWidth();
        double strHeight = bounds.getHeight();
        double strPosition = (w - strWidth) / 2;

        double lineY = getPrefHeight() - 35 - strHeight / 2;

        graphicsContext.strokeLine(0, lineY, (w - strWidth) / 2 - 10, lineY);
        double[] arrowX = {0.0, 10.0, 10.0};
        double[] arrowY = {lineY, lineY + 3, lineY - 3};
        graphicsContext.fillPolygon(arrowX, arrowY, arrowX.length);

        graphicsContext.strokeLine((w + strWidth) / 2 + 10, lineY, w, lineY);
        arrowX = new double[]{w, w - 10, w - 10};
        graphicsContext.fillPolygon(arrowX, arrowY, arrowX.length);

        graphicsContext.fillText(rangeString, strPosition, getPrefHeight() - 35);

    }

    private void drawTicks(GraphicsContext graphicsContext, InsertionMarker i) {

        double w = getPrefWidth();
        if (w < 200) {
            return;
        }

        graphicsContext.setFont(tickFont);

        // location unit for whole genome is kilobase
        boolean scaleInKB = frame.getChrName().equals(Globals.CHR_ALL);

        double range = w * frame.getScale();
        TickSpacing ts = findSpacing(range, scaleInKB);
        double spacing = ts.getMajorTick();

        // Find starting point closest to the current origin
        int nTick = (int) (frame.getOrigin() / spacing) - 1;
        double l = (int) (nTick * spacing);
        double x = frame.getScreenPosition(l - 1 + 0.5);    // 0 vs 1 based coordinates, then center over base
        //int strEnd = Integer.MIN_VALUE;

        Text sizer = new Text("");
        sizer.setFont(graphicsContext.getFont());
        while (x < getPrefWidth()) {
            l = (int) (nTick * spacing);
            x = frame.getScreenPosition(l - 1 + 0.5);
            String chrPosition = formatNumber((double) l / ts.getUnitMultiplier()) +
                    " " + ts.getMajorUnit();
            double strWidth = FontMetrics.getTextWidthInFont(chrPosition, sizer);
            double strPosition = x - strWidth / 2;
            //if (strPosition > strEnd) {

            final double height = getPrefHeight();
            if (nTick % 2 == 0) {
                graphicsContext.fillText(chrPosition, strPosition, height - 15);
            }

            graphicsContext.strokeLine(x, height - 10, x, height - 2);
            nTick++;
        }
    }

    private void drawChromosomeTicks(GraphicsContext graphicsContext) {

        // TODO -- remove hardcoded value
        int locationUnit = 1000;

        graphicsContext.setFont(chrFont);

        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        if (genome == null) {
            log.warn("No genome found");
            PreferencesManager.getPreferences().remove(DEFAULT_GENOME);
            return;
        }

        boolean even = true;
        long offset = 0;
        chromosomeRects.clear();
        List<String> chrNames = genome.getLongChromosomeNames();
        if (chrNames == null) {
            log.info("No chromosomes found for genome: " + PreferencesManager.getPreferences().getDefaultGenome());
            PreferencesManager.getPreferences().remove(DEFAULT_GENOME);
        }
        if (chrNames.size() > 500) {
            return;
        }

        Text sizer = new Text("");
        sizer.setFont(graphicsContext.getFont());
        for (String chrName : chrNames) {
            Chromosome c = genome.getChromosome(chrName);
            if (c == null) {
                log.info("Chromosome '" + chrName + "' not found");
                continue;
            }
            int chrLength = c.getLength();

            double scale = frame.getScale();
            int gStart = genome.getGenomeCoordinate(chrName, 0);
            double x = gStart / scale;
            double dw = chrLength / (locationUnit * scale);


            graphicsContext.strokeLine(x, getPrefHeight() - 10, x, getPrefHeight() - 2);

            // Don't label chromosome if its width is < 5 pixels
            if (dw > 5) {
                double center = x + dw / 2;

                String displayName = null;
                if (chrName.startsWith("gi|")) {
                    displayName = Genome.getNCBIName(chrName);
                } else {
                    displayName = chrName.replace("chr", "");
                }
                double strWidth = FontMetrics.getTextWidthInFont(displayName, sizer);
                double strPosition = center - strWidth / 2;


                double y = (even ? getPrefHeight() - 35 : getPrefHeight() - 25);
                graphicsContext.fillText(displayName, strPosition, y);
                Rectangle2D clickRect = new Rectangle2D(strPosition, y - 15, strWidth, 15);
                String tooltipText = "Jump to chromosome: " + chrName;
                chromosomeRects.add(new ClickLink(clickRect, chrName, tooltipText));

                even = !even;

            }

            offset += chrLength;
        }
    }

    public static String formatNumber(double position) {
        return DECIMAL_FORMAT.format((int) position);
    }

    public static TickSpacing findSpacing(double maxValue, boolean scaleInKB) {

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
        // TODO: for now not dealing with Mouse events, cursors, DnD, tooltips, etc

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

        Rectangle2D region;
        String value;
        String tooltipText;

        ClickLink(Rectangle2D region, String value, String tooltipText) {
            this.region = region;
            this.value = value;
            this.tooltipText = tooltipText;
        }
    }


    static class MouseRect {
        Rectangle2D bounds;
        String text;

        MouseRect(Rectangle2D bounds, String text) {
            this.bounds = bounds;
            this.text = text;
        }

        boolean containsPoint(Point2D p) {
            return bounds.contains(p);
        }

        String getText() {
            return text;
        }

        public double width() {
            return bounds.getWidth();
        }
    }


}
