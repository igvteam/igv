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
package org.broad.igv.sam;

import htsjdk.samtools.util.SequenceUtil;
import org.broad.igv.Globals;
import org.broad.igv.data.CoverageDataSource;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.renderer.*;
import org.broad.igv.sam.mods.BaseModificationCounts;
import org.broad.igv.sam.mods.BaseModificationUtils;
import org.broad.igv.tdf.TDFDataSource;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.track.*;
import org.broad.igv.ui.DataRangeDialog;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.StringUtils;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import static org.broad.igv.prefs.Constants.*;

/**
 * @author jrobinso
 */

public class CoverageTrack extends AbstractTrack implements ScalableTrack {

    private static Logger log = LogManager.getLogger(CoverageTrack.class);

    public static final int TEN_MB = 10000000;
    static DecimalFormat locationFormatter = new DecimalFormat();

    char[] nucleotides = {'a', 'c', 'g', 't', 'n'};

    public static final boolean DEFAULT_AUTOSCALE = true;

    private Color _color = null;   // Explicit color setting
    private float snpThreshold;
    private AlignmentTrack alignmentTrack;
    private AlignmentDataManager dataManager;
    private CoverageDataSource dataSource;
    private DataRenderer dataSourceRenderer;
    private IntervalRenderer intervalRenderer;
    private IGVPreferences prefs;
    private Genome genome;
    private boolean removed = false;
    IGV igv;

    ColorScale baseModificationColorScale;


    /**
     * Whether to autoscale across all ReferenceFrames
     * Default is true because we usually do, SashimiPlot does not
     */
    private boolean globalAutoScale = true;

    public CoverageTrack() {
    }

    /**
     * Copy constructor.  Used for Sashimi plot.
     *
     * @param track
     */
    public CoverageTrack(CoverageTrack track) {
        this(track.getResourceLocator(), track.getName(), track.alignmentTrack, track.genome);
        if (track.dataManager != null) this.setDataManager(track.dataManager);
        if (track.dataSource != null) this.setDataSource(track.dataSource);
        this.snpThreshold = track.snpThreshold;
        this.prefs = track.prefs;
        this.igv = IGV.hasInstance() ? IGV.getInstance() : null;
    }

    public CoverageTrack(ResourceLocator locator, String name, AlignmentTrack alignmentTrack, Genome genome) {
        super(locator, locator.getPath() + "_coverage", name);
        super.setDataRange(new DataRange(0, 0, 60));
        this.alignmentTrack = alignmentTrack;
        this.genome = genome;
        intervalRenderer = new IntervalRenderer();
        prefs = PreferencesManager.getPreferences();
        snpThreshold = prefs.getAsFloat(SAM_ALLELE_THRESHOLD);
        autoScale = DEFAULT_AUTOSCALE;
        this.igv = IGV.hasInstance() ? IGV.getInstance() : null;
    }

    @Override
    public Color getColor() {
        return _color == null ?
                ColorUtilities.slightlyDarker(alignmentTrack.getColor()) :
                _color;
    }

    @Override
    public void setColor(Color color) {
        _color = color;
    }

    @Override
    public String getSample() {
        if (sampleId != null) {
            return sampleId;    // Explicitly set sample ID (e.g. from server load XML)
        }
        return alignmentTrack == null ? null : alignmentTrack.getSample();
    }

    @Override
    public boolean isNumeric() {
        return true;
    }

    public void setDataManager(AlignmentDataManager dataManager) {
        this.dataManager = dataManager;
        this.dataManager.subscribe(this);
    }

    public void setDataSource(CoverageDataSource dataSource) {
        this.dataSource = dataSource;
        dataSourceRenderer = new BarChartRenderer();
        setDataRange(new DataRange(0, 0, 1.5f * (float) dataSource.getDataMax()));

    }

    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {
        if (frame.getChrName().equals(Globals.CHR_ALL) || frame.getScale() > dataManager.getMinVisibleScale()) {
            return true;   // Nothing to paint
        } else {
            return dataManager.isLoaded(frame);
        }
    }

    @Override
    public void load(ReferenceFrame referenceFrame) {
        dataManager.load(referenceFrame, alignmentTrack.renderOptions, true);
    }


    public void setSnpThreshold(float snpThreshold) {
        this.snpThreshold = snpThreshold;
    }

    public float getSnpThreshold() {
        return snpThreshold;
    }

    public boolean isRemoved() {
        return removed;
    }

    @Override
    public boolean isVisible() {
        return super.isVisible() && !removed;
    }

    @Override
    public void unload() {
        super.unload();
        removed = true;
        if (dataManager != null) {
            dataManager.unsubscribe(this);
        }
        setVisible(false);
    }

    public void render(RenderContext context, Rectangle rect) {

        int viewWindowSize = context.getReferenceFrame().getCurrentRange().getLength();
        if (viewWindowSize > dataManager.getVisibilityWindow() && dataSource == null) {
            Rectangle visibleRect = context.getVisibleRect().intersection(rect);
            Graphics2D g = context.getGraphic2DForColor(Color.gray);
            GraphicUtils.drawCenteredText("Zoom in to see coverage.", visibleRect, g);
            return;
        }


        drawData(context, rect);

        drawBorder(context, rect);

        if (dataSourceRenderer != null) {
            dataSourceRenderer.renderBorder(this, context, rect);
            dataSourceRenderer.renderAxis(this, context, rect);
        } else {
            DataRenderer.drawScale(this.getDataRange(), context, rect);
        }
    }


    public void drawData(RenderContext context, Rectangle rect) {

        int viewWindowSize = context.getReferenceFrame().getCurrentRange().getLength();
        if (viewWindowSize <= dataManager.getVisibilityWindow() && !context.getChr().equals(Globals.CHR_ALL)) {
            //Show coverage calculated from intervals if zoomed in enough
            AlignmentInterval interval = null;
            if (dataManager != null) {
                interval = dataManager.getLoadedInterval(context.getReferenceFrame(), true);
            }
            if (interval != null) {
                intervalRenderer.paint(context, rect, interval.getCounts());
                return;
            }
        }

        //Not rendered yet.  Use precomputed scores, if available
        List<LocusScore> scores = getInViewScores(context.getReferenceFrame());
        if (scores != null) {
            dataSourceRenderer.renderScores(this, scores, context, rect);
        }

    }


    private List<LocusScore> getInViewScores(ReferenceFrame frame) {

        List<LocusScore> inViewScores = null;

        if (dataSource != null) {
            String chr = frame.getChrName();
            int start = (int) frame.getOrigin();
            int end = (int) frame.getEnd();
            int zoom = frame.getZoom();
            inViewScores = dataSource.getSummaryScoresForRange(chr, start, end, zoom);

            // Trim
            // Trim scores
            int startIdx = Math.max(0, FeatureUtils.getIndexBefore(start, inViewScores));
            int endIdx = inViewScores.size() - 1;   // Starting guess
            int tmp = Math.max(0, FeatureUtils.getIndexBefore(end, inViewScores));
            for (int i = tmp; i < inViewScores.size(); i++) {
                if (inViewScores.get(i).getStart() > end) {
                    endIdx = i - 1;
                    break;
                }
            }
            endIdx = Math.max(startIdx + 1, endIdx);

            if (inViewScores.size() > 1) {
                return startIdx == 0 && endIdx == inViewScores.size() - 1 ?
                        inViewScores :
                        inViewScores.subList(startIdx, endIdx);
            } else {
                return inViewScores;
            }
        }
        return inViewScores;
    }


    @Override
    public Range getInViewRange(ReferenceFrame frame) {

        int viewWindowSize = frame.getCurrentRange().getLength();
        if (dataManager == null || viewWindowSize > dataManager.getVisibilityWindow()) {
            List<LocusScore> scores = getInViewScores(frame);
            if (scores != null && scores.size() > 0) {
                float min = scores.get(0).getScore();
                float max = min;
                for (int i = 1; i < scores.size(); i++) {
                    LocusScore score = scores.get(i);
                    float value = score.getScore();
                    min = Math.min(value, min);
                    max = Math.max(value, max);
                }
                return new Range(min, max);
            } else {
                return null;
            }

        } else {
            AlignmentInterval interval = dataManager.getLoadedInterval(frame);
            if (interval == null) return null;

            int origin = (int) frame.getOrigin();
            int end = (int) frame.getEnd() + 1;
            int intervalMax = interval.getMaxCount(origin, end);

            return new Range(0, Math.max(10, intervalMax));
        }
    }


    /**
     * Draw border and scale
     *
     * @param context
     * @param rect
     */
    private void drawBorder(RenderContext context, Rectangle rect) {
        context.getGraphic2DForColor(Color.gray).drawLine(
                rect.x, rect.y + rect.height,
                rect.x + rect.width, rect.y + rect.height);
    }

    public void drawScale(RenderContext context, Rectangle rect) {
        DataRenderer.drawScale(getDataRange(), context, rect);
    }

    public boolean isLogNormalized() {
        return false;
    }

    public String getValueStringAt(String chr, double position, int mouseX, int mouseY, ReferenceFrame frame) {

        float maxRange = PreferencesManager.getPreferences().getAsFloat(SAM_MAX_VISIBLE_RANGE);
        float minVisibleScale = (maxRange * 1000) / 700;

        StringBuffer buf = new StringBuffer();

        if (!chr.equals("All")) {
            String posString = chr + ":" + locationFormatter.format(Math.floor(position + 1));
            buf.append(posString + "<br>");
            buf.append("<hr>");
        }

        if (frame.getScale() < minVisibleScale) {
            AlignmentInterval interval = dataManager.getLoadedInterval(frame);
            if (interval != null && interval.contains(chr, (int) position, (int) position)) {
                AlignmentCounts counts = interval.getCounts();
                if (counts != null) {
                    buf.append(counts.getValueStringAt((int) position));
                    final AlignmentTrack.ColorOption colorOption = alignmentTrack.renderOptions.getColorOption();
                    if ((colorOption == AlignmentTrack.ColorOption.BASE_MODIFICATION ||
                            colorOption == AlignmentTrack.ColorOption.BASE_MODIFICATION_5MC) &&
                            counts.getModifiedBaseCounts() != null) {
                        buf.append("<hr>");
                        buf.append(counts.getModifiedBaseCounts().getValueString((int) position, colorOption));
                    }
                }
            }
        } else {
            buf.append(getPrecomputedValueString(chr, position, frame));
        }
        return buf.toString();
    }

    private String getPrecomputedValueString(String chr, double position, ReferenceFrame frame) {

        if (dataSource == null) {
            return "";
        }
        int zoom = Math.max(0, frame.getZoom());
        List<LocusScore> scores = dataSource.getSummaryScoresForRange(chr, (int) position - 10, (int) position + 10, zoom);

        // give a 2 pixel window, otherwise very narrow features will be missed.
        double bpPerPixel = frame.getScale();
        double minWidth = 2 * bpPerPixel;    /* * */

        if (scores == null) {
            return "";
        } else {
            LocusScore score = (LocusScore) FeatureUtils.getFeatureAt(position, 0, scores);
            return score == null ? "" : "Mean count: " + score.getScore();
        }
    }

    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, String frameName) {
        return 0;
    }

    public void rescale(ReferenceFrame iframe) {
        List<ReferenceFrame> frameList = new ArrayList<ReferenceFrame>();
        if (iframe != null) frameList.add(iframe);
        if (globalAutoScale) {
            frameList.addAll(FrameManager.getFrames());
        }

        if (autoScale && dataManager != null) {

            int max = 10;
            for (ReferenceFrame frame : frameList) {
                AlignmentInterval interval = dataManager.getLoadedInterval(frame);
                if (interval == null) continue;

                int origin = (int) frame.getOrigin();
                int end = (int) frame.getEnd() + 1;

                int intervalMax = interval.getMaxCount(origin, end);
                max = intervalMax > max ? intervalMax : max;
            }

            DataRange newRange = new DataRange(0, max);
            newRange.setType(getDataRange().getType());
            super.setDataRange(newRange);

        }
    }

    /**
     * Class to render coverage track, including mismatches.
     * <p/>
     * NOTE:  This class has been extensively optimized with the aid of a profiler,  attempts to "clean up" this code
     * should be done with frequent profiling, or it will likely have detrimental performance impacts.
     */

    /**
     *
     */
    class IntervalRenderer {

        /**
         * @param context
         * @param rect            - the track rectangle
         * @param alignmentCounts
         */
        private void paint(final RenderContext context, final Rectangle rect, final AlignmentCounts alignmentCounts) {

            Color color = getColor();
            Graphics2D graphics = context.getGraphic2DForColor(color);
            final AlignmentTrack.ColorOption colorOption = alignmentTrack.renderOptions.getColorOption();
            boolean bisulfiteMode = colorOption == AlignmentTrack.ColorOption.BISULFITE;
            boolean baseModMode = colorOption == AlignmentTrack.ColorOption.BASE_MODIFICATION ||
                    colorOption == AlignmentTrack.ColorOption.BASE_MODIFICATION_5MC;

            final int intervalEnd = alignmentCounts.getEnd();
            final int intervalStart = alignmentCounts.getStart();

            BisulfiteCounts bisulfiteCounts = null;
            byte[] refBases = null;
            if ((intervalEnd - intervalStart) < TEN_MB && alignmentCounts.hasBaseCounts()) {
                refBases = genome.getSequence(context.getChr(), intervalStart, intervalEnd);
                bisulfiteCounts = alignmentCounts.getBisulfiteCounts();
            }

            DataRange range = getDataRange();
            double maxRange = range.isLog() ? Math.log10(range.getMaximum() + 1) : range.getMaximum();
            final double origin = context.getOrigin();
            final double scale = context.getScale();
            int start = alignmentCounts.getStart();
            int step = alignmentCounts.getBucketSize();
            int nPoints = alignmentCounts.getNumberOfPoints();
            boolean isSparse = alignmentCounts instanceof SparseAlignmentCounts;

            // First pass -- draw gray coverage bars
            for (int idx = 0; idx < nPoints; idx++) {

                int pos = isSparse ? ((SparseAlignmentCounts) alignmentCounts).getPosition(idx) : start + idx * step;
                int pX = (int) (rect.x + (pos - origin) / scale);
                double endX = rect.x + (pos + step - origin) / scale;
                int dX = (int) ((endX - pX) < 1 ? 1 : (endX - pX) > 3 ? endX - pX - 1 : endX - pX);

                if (pX > rect.x + rect.width) {
                    break; // We're done,  data is position sorted so we're beyond the right-side of the view
                } else if (endX < rect.x) {
                    continue;
                }

                int totalCount = alignmentCounts.getTotalCount(pos);
                double tmp = range.isLog() ? Math.log10(totalCount + 1) / maxRange : totalCount / maxRange;
                int barHeight = (int) Math.min(tmp * rect.height, rect.height - 1);
                if (barHeight > 0) {
                    int bottomY = rect.y + rect.height;
                    int topY = bottomY - barHeight;
                    graphics.fillRect(pX, topY, dX, barHeight);
                }
            }

            // Second pass -- potentially overlay mismatches
            for (int idx = 0; idx < nPoints; idx++) {

                int pos = isSparse ? ((SparseAlignmentCounts) alignmentCounts).getPosition(idx) : start + idx * step;
                int pX = (int) (rect.x + (pos - origin) / scale);
                double endX = rect.x + (pos + step - origin) / scale;
                int dX = (int) ((endX - pX) < 1 ? 1 : (endX - pX) > 3 ? endX - pX - 1 : endX - pX);

                if (pX > rect.x + rect.width) {
                    break; // We're done,  data is position sorted so we're beyond the right-side of the view
                } else if (endX < rect.x) {
                    continue;
                }

                int totalCount = alignmentCounts.getTotalCount(pos);
                double tmp = range.isLog() ? Math.log10(totalCount + 1) / maxRange : totalCount / maxRange;
                int barHeight = (int) Math.min(tmp * rect.height, rect.height - 1);
                if (barHeight > 0) {
                    int bottomY = rect.y + rect.height;

                    // Potentially color mismatch
                    if (bisulfiteMode) {
                        BisulfiteCounts.Count bc = bisulfiteCounts != null ? bisulfiteCounts.getCount(pos) : null;
                        if (bc != null && (bc.methylatedCount + bc.unmethylatedCount) > 0) {
                            drawBarBisulfite(context, pX, bottomY, dX, barHeight, totalCount, bc);
                        }
                    } else if (baseModMode) {
                        drawModifiedBaseBar(context, pX, bottomY, dX, barHeight, pos, alignmentCounts);
                    } else {
                        if (refBases != null) {
                            int refIdx = pos - intervalStart;
                            if (refIdx >= 0 && refIdx < refBases.length) {
                                byte ref = refBases[refIdx];
                                if (alignmentCounts.isConsensusMismatch(pos, ref, context.getChr(), snpThreshold)) {
                                    drawAllelFreqBar(context, pX, bottomY, dX, barHeight, pos, totalCount, alignmentCounts);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * Draw a colored bar to represent a mismatch to the reference.   The height is proportional to the % of
     * reads with respect to the total.  If "showAllSnps == true"  the bar is shaded by avg read quality.
     *
     * @param context
     * @param pos
     * @param barHeight       -- total height in pixels of the coverage bar
     * @param pBottom
     * @param pX
     * @param dX
     * @param alignmentCounts
     * @return
     */

    void drawAllelFreqBar(RenderContext context,
                          int pX,
                          int pBottom,
                          int dX,
                          int barHeight,
                          int pos,
                          double totalCount,
                          AlignmentCounts alignmentCounts) {

        for (char nucleotide : nucleotides) {
            int count = alignmentCounts.getCount(pos, (byte) nucleotide);
            if (count > 0) {
                double f = count / totalCount;
                int alleleHeight = (int) (f * barHeight);
                int baseY = pBottom - alleleHeight;
                Color c = SequenceRenderer.nucleotideColors.get(nucleotide);
                Graphics2D tGraphics = context.getGraphic2DForColor(c);
                tGraphics.fillRect(pX, baseY, dX, alleleHeight);
                pBottom = baseY;
            }
        }
    }

    int drawModifiedBaseBar(RenderContext context,
                            int pX,
                            int pBottom,
                            int dX,
                            int barHeight,
                            int pos,
                            AlignmentCounts alignmentCounts) {

        BaseModificationCounts baseCounts = alignmentCounts.getModifiedBaseCounts();
        AlignmentTrack.ColorOption colorOption = alignmentTrack.renderOptions.getColorOption();
        if (baseCounts != null) {

            Graphics2D graphics = context.getGraphics();

            for (BaseModificationCounts.Key key : baseCounts.getAllModifications()) {

                int count = baseCounts.getCount(pos, key, colorOption);

                if (barHeight > 0 && count > 0) {

                    byte base = (byte) key.getBase();
                    byte complement = SequenceUtil.complement(base);
                    char modStrand = key.getStrand();
                    String modification = key.getModification();

                    int cCount = modStrand == '+' ?
                            alignmentCounts.getPosCount(pos, base) + alignmentCounts.getNegCount(pos, complement) :
                            alignmentCounts.getPosCount(pos, complement) + alignmentCounts.getNegCount(pos, base);

                    int calledBarHeight = (int) ((((float) count) / cCount) * barHeight);
                    Color noModColor = BaseModificationUtils.getModColor(modification, (byte) 0, colorOption);
                    Color modColor = BaseModificationUtils.getModColor(modification, (byte) 255, colorOption);

                    float averageLikelihood = (float) (baseCounts.getLikelhoodSum(pos, key)) / (count * 255);
                    int modHeight = (int) (averageLikelihood * calledBarHeight);

                    if (colorOption == AlignmentTrack.ColorOption.BASE_MODIFICATION_5MC) {

                        int noModHeight = calledBarHeight - modHeight;
                        int baseY = pBottom - noModHeight;

                        graphics.setColor(noModColor);
                        graphics.fillRect(pX, baseY, dX, noModHeight);

                        if(modHeight > 0) {
                            baseY -= modHeight;
                            graphics.setColor(modColor);
                            graphics.fillRect(pX, baseY, dX, modHeight);
                        }

                    } else {
                        // Generic modification
                        float threshold = PreferencesManager.getPreferences().getAsFloat("SAM.BASEMOD_THRESHOLD");
                        if(averageLikelihood > threshold && modHeight > 0) {
                            int baseY = pBottom - modHeight;
                            graphics.setColor(modColor);
                            graphics.fillRect(pX, baseY, dX, modHeight);
                            pBottom = baseY;
                        }
                    }

                }
            }
        }
        return pX + dX;
    }


    void drawBarBisulfite(RenderContext context,
                          int pX0,
                          int pBottom,
                          int dX,
                          int barHeight,
                          double totalCount,
                          BisulfiteCounts.Count count) {

        // If bisulfite mode, we expand the rectangle to make it more visible.  This code is copied from AlignmentRenderer
        int pX = pX0;
        if (dX < 3) {
            int expansion = dX;
            pX -= expansion;
            dX += (2 * expansion);
        }


        double nMethylated = count.methylatedCount;
        double unMethylated = count.unmethylatedCount;
        Color c = Color.red;
        Graphics2D tGraphics = context.getGraphic2DForColor(c);

        //Not all reads at a position are informative,  color by % of informative reads
        // double totalInformative = count.methylatedCount + count.unmethylatedCount;
        // double mult = totalCount / totalInformative;
        // nMethylated *= mult;
        // unMethylated *= mult;

        double f = nMethylated / totalCount;
        int height = (int) (f * barHeight);

        int baseY = pBottom - height;
        if (height > 0) {
            tGraphics.fillRect(pX, baseY, dX, height);
        }
        pBottom = baseY;

        c = Color.blue;
        tGraphics = context.getGraphic2DForColor(c);

        f = unMethylated / totalCount;
        height = (int) (f * barHeight);
        if (height > 0) {
            baseY = pBottom - height;
            tGraphics.fillRect(pX, baseY, dX, height);
        }
    }


    /**
     * Strand-specific
     *
     * @param context
     * @param pos
     * @param rect
     * @param maxCount
     * @param pY
     * @param pX
     * @param dX
     * @param isPositive
     * @param interval
     * @return
     */
    void drawStrandBar(RenderContext context,
                       int pos,
                       Rectangle rect,
                       double maxCount,
                       int pY,
                       int pX,
                       int dX,
                       boolean isPositive,
                       AlignmentCounts interval) {


        for (char nucleotide : nucleotides) {

            Color c = SequenceRenderer.nucleotideColors.get(nucleotide);
            Graphics2D tGraphics = context.getGraphic2DForColor(c);

            int count = isPositive ? interval.getPosCount(pos, (byte) nucleotide) :
                    interval.getNegCount(pos, (byte) nucleotide);

            int height = (int) Math.round(count * rect.getHeight() / maxCount);
            height = isPositive ? Math.min(pY - rect.y, height) :
                    Math.min(rect.y + rect.height - pY, height);
            int baseY = (int) (isPositive ? (pY - height) : pY);

            if (height > 0) {
                tGraphics.fillRect(pX, baseY, dX, height);
            }
            pY = isPositive ? baseY : baseY + height;
        }
    }


    static float[] colorComps = new float[3];

    private Color getShadedColor(int qual, Color backgroundColor, Color color) {
        float alpha = 0;
        int minQ = prefs.getAsInt(SAM_BASE_QUALITY_MIN);
        ColorUtilities.getRGBColorComponents(color);
        if (qual < minQ) {
            alpha = 0.1f;
        } else {
            int maxQ = prefs.getAsInt(SAM_BASE_QUALITY_MAX);
            alpha = Math.max(0.1f, Math.min(1.0f, 0.1f + 0.9f * (qual - minQ) / (maxQ - minQ)));
        }
        // Round alpha to nearest 0.1, for effeciency;
        alpha = ((int) (alpha * 10 + 0.5f)) / 10.0f;

        if (alpha >= 1) {
            return color;
        } else {
            return ColorUtilities.getCompositeColor(backgroundColor, color, alpha);
        }
    }

    /**
     * Override to return a specialized popup menu
     *
     * @return
     */
    @Override
    public IGVPopupMenu getPopupMenu(TrackClickEvent te) {

        Collection<Track> tmp = new ArrayList<Track>();
        tmp.add(this);

        IGVPopupMenu popupMenu = TrackMenuUtils.getPopupMenu(tmp, getName(), te);

        popupMenu.addSeparator();
        this.addSnpTresholdItem(popupMenu);

        popupMenu.addSeparator();
        addLoadCoverageDataItem(popupMenu);

        popupMenu.addSeparator();
        addCopyDetailsItem(popupMenu, te);

        if (alignmentTrack != null) {
            popupMenu.addSeparator();
            addShowItems(popupMenu);
        }

        return popupMenu;
    }

    private void addCopyDetailsItem(IGVPopupMenu popupMenu, TrackClickEvent te) {
        JMenuItem copyDetails = new JMenuItem("Copy Details to Clipboard");
        copyDetails.setEnabled(false);
        if (te.getFrame() != null) {
            final String details = getValueStringAt(te.getFrame().getChrName(), te.getChromosomePosition(), te.getMouseEvent().getX(), te.getMouseEvent().getY(), te.getFrame());
            copyDetails.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    if (details != null) {
                        String deets = details.replace("<br>", System.getProperty("line.separator"));
                        StringUtils.copyTextToClipboard(deets);
                    }
                }
            });
            copyDetails.setEnabled(details != null);
        }
        popupMenu.add(copyDetails);
    }


    public static JMenuItem addDataRangeItem(final Frame parentFrame, JPopupMenu menu, final Collection<? extends Track> selectedTracks) {
        JMenuItem maxValItem = new JMenuItem("Set Data Range");

        maxValItem.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                if (selectedTracks.size() > 0) {

                    DataRange prevAxisDefinition = selectedTracks.iterator().next().getDataRange();
                    DataRangeDialog dlg = new DataRangeDialog(parentFrame, prevAxisDefinition);
                    dlg.setHideMid(true);
                    dlg.setVisible(true);
                    if (!dlg.isCanceled()) {
                        float min = Math.min(dlg.getMin(), dlg.getMax());
                        float max = Math.max(dlg.getMin(), dlg.getMax());
                        float mid = dlg.getBase();
                        if (mid < min) mid = min;
                        else if (mid > max) mid = max;
                        DataRange dataRange = new DataRange(min, mid, max);
                        dataRange.setType(dlg.getDataRangeType());

                        for (Track track : selectedTracks) {
                            track.setDataRange(dataRange);
                            track.setAutoScale(false);
                        }
                        parentFrame.repaint();
                    }
                }

            }
        });
        if (menu != null) menu.add(maxValItem);

        return maxValItem;
    }

    public JMenuItem addSnpTresholdItem(JPopupMenu menu) {
        JMenuItem maxValItem = new JMenuItem("Set allele frequency threshold...");
        maxValItem.addActionListener(e -> {
            String value = JOptionPane.showInputDialog("Allele frequency threshold: ", Float.valueOf(snpThreshold));
            if (value == null) {
                return;
            }
            try {
                float tmp = Float.parseFloat(value);
                snpThreshold = tmp;
                repaint();
            } catch (Exception exc) {
                //log
            }
        });
        menu.add(maxValItem);
        return maxValItem;
    }

    public void addShowItems(JPopupMenu menu) {

        final SpliceJunctionTrack spliceJunctionTrack = alignmentTrack.getSpliceJunctionTrack();

        final JMenuItem item = new JCheckBoxMenuItem("Show Coverage Track");
        item.setSelected(true);
        item.addActionListener(e -> {
            CoverageTrack.this.setVisible(item.isSelected());
            IGV.getInstance().repaint(Arrays.asList(CoverageTrack.this));
        });
        //If this is the only track visible, disable option to hide it.
        if (!(alignmentTrack.isVisible() ||
                (spliceJunctionTrack != null && spliceJunctionTrack.isVisible()))) {
            item.setEnabled(false);
        }
        menu.add(item);

        if (spliceJunctionTrack != null) {
            final JMenuItem junctionItem = new JCheckBoxMenuItem("Show Splice Junction Track");
            junctionItem.setSelected(spliceJunctionTrack.isVisible());
            junctionItem.setEnabled(!spliceJunctionTrack.isRemoved());
            junctionItem.addActionListener(e -> {
                spliceJunctionTrack.setVisible(junctionItem.isSelected());
                IGV.getInstance().repaint(Arrays.asList(spliceJunctionTrack));
            });
            menu.add(junctionItem);
        }

        if (alignmentTrack != null) {
            final JMenuItem alignmentItem = new JCheckBoxMenuItem("Show Alignment Track");
            alignmentItem.setSelected(alignmentTrack.isVisible());
            alignmentItem.setEnabled(!alignmentTrack.isRemoved());
            alignmentItem.addActionListener(e -> {
                alignmentTrack.setVisible(alignmentItem.isSelected());
                IGV.getInstance().repaint(Arrays.asList(alignmentTrack));
            });
            menu.add(alignmentItem);
        }
    }


    public void addLoadCoverageDataItem(JPopupMenu menu) {
        // Change track height by attribute
        final JMenuItem item = new JCheckBoxMenuItem("Load pre-computed coverage data...");
        item.addActionListener(e -> {
            final IGVPreferences prefs = PreferencesManager.getPreferences();
            File initDirectory = prefs.getLastTrackDirectory();
            File file = FileDialogUtils.chooseFile("Select coverage file", initDirectory, FileDialog.LOAD);
            if (file != null) {
                prefs.setLastTrackDirectory(file.getParentFile());
                String path = file.getAbsolutePath();
                if (path.endsWith(".tdf") || path.endsWith(".tdf")) {
                    TDFReader reader = TDFReader.getReader(file.getAbsolutePath());
                    TDFDataSource ds = new TDFDataSource(reader, 0, getName() + " coverage", genome);
                    setDataSource(ds);
                    repaint();
                } else {
                    if (igv != null) MessageUtils.showMessage("Coverage data must be in .tdf format");
                }
            }
        });

        item.setEnabled(dataSource == null);
        menu.add(item);

    }


    public void setGlobalAutoScale(boolean globalAutoScale) {
        this.globalAutoScale = globalAutoScale;
    }


    @Override
    public void marshalXML(Document document, Element element) {

        super.marshalXML(document, element);
        element.setAttribute("snpThreshold", String.valueOf(snpThreshold));
        if (_color != null) {
            element.setAttribute("_color", ColorUtilities.colorToString(_color));
        }

    }

    @Override
    public void unmarshalXML(Element element, Integer version) {
        super.unmarshalXML(element, version);
        if (element.hasAttribute("snpThreshold")) {
            snpThreshold = Float.parseFloat(element.getAttribute("snpThreshold"));
        }
        if (element.hasAttribute("_color")) {
            _color = ColorUtilities.stringToColor(element.getAttribute("_color"));
        }

    }

}