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

import com.google.common.eventbus.Subscribe;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.data.CoverageDataSource;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.goby.GobyCountArchiveDataSource;
import org.broad.igv.renderer.BarChartRenderer;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.DataRenderer;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.session.IGVSessionReader;
import org.broad.igv.session.SubtlyImportant;
import org.broad.igv.tdf.TDFDataSource;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.track.*;
import org.broad.igv.ui.DataRangeDialog;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.event.AlignmentTrackEvent;
import org.broad.igv.ui.event.DataLoadedEvent;
import org.broad.igv.ui.event.ViewChange;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.StringUtils;

import javax.swing.*;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlType;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

/**
 * @author jrobinso
 */
@XmlType(factoryMethod = "getNextTrack")

public class CoverageTrack extends AbstractTrack implements ScalableTrack {

    private static Logger log = Logger.getLogger(CoverageTrack.class);
    public static final int TEN_MB = 10000000;

    char[] nucleotides = {'a', 'c', 'g', 't', 'n'};
    public static Color lightBlue = new Color(0, 0, 150);
    private static Color coverageGrey = new Color(175, 175, 175);
    public static final Color negStrandColor = new Color(140, 140, 160);
    public static final Color posStrandColor = new Color(160, 140, 140);

    public static final boolean DEFAULT_AUTOSCALE = true;
    public static final boolean DEFAULT_SHOW_REFERENCE = false;

    // User settable state -- these attributes should be stored in the session file
    @XmlAttribute
    private float snpThreshold;
    //TODO This appears to not be used anywhere, remove?
    @XmlAttribute
    private boolean showReference;

    private AlignmentTrack alignmentTrack;
    private AlignmentDataManager dataManager;
    private CoverageDataSource dataSource;
    private DataRenderer dataSourceRenderer;
    private IntervalRenderer intervalRenderer;
    private PreferenceManager prefs;
    private JMenuItem dataRangeItem;
    private Genome genome;
    private boolean removed = false;

    /**
     * Whether to autoscale across all ReferenceFrames
     * Default is true because we usually do, SashimiPlot does not
     */
    private boolean globalAutoScale = true;

    private AlignmentTrack.RenderOptions renderOptions = null;

    @SubtlyImportant
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
    }

    public CoverageTrack(ResourceLocator locator, String name, AlignmentTrack alignmentTrack, Genome genome) {
        super(locator, locator.getPath() + "_coverage", name);
        super.setDataRange(new DataRange(0, 0, 60));
        this.alignmentTrack = alignmentTrack;
        this.genome = genome;
        intervalRenderer = new IntervalRenderer();
        setMaximumHeight(40);

        setColor(coverageGrey);

        prefs = PreferenceManager.getInstance();
        snpThreshold = prefs.getAsFloat(PreferenceManager.SAM_ALLELE_THRESHOLD);
        autoScale = DEFAULT_AUTOSCALE;
        showReference = DEFAULT_SHOW_REFERENCE;
        //TODO  logScale = prefs.


    }

    public void setDataManager(AlignmentDataManager dataManager) {
        this.dataManager = dataManager;
        this.dataManager.getEventBus().register(this);
    }

    public void setDataSource(CoverageDataSource dataSource) {
        this.dataSource = dataSource;
        dataSourceRenderer = new BarChartRenderer();
        setDataRange(new DataRange(0, 0, 1.5f * (float) dataSource.getDataMax()));

    }

    public void setSnpThreshold(float snpThreshold) {
        this.snpThreshold = snpThreshold;
    }

    public float getSnpThreshold() {
        return snpThreshold;
    }

    public boolean isShowReference() {
        return showReference;
    }


    public void setRenderOptions(AlignmentTrack.RenderOptions renderOptions) {
        this.renderOptions = renderOptions;
    }

    AlignmentTrack.RenderOptions getRenderOptions() {
        return this.renderOptions;
    }

    public boolean isRemoved() {
        return removed;
    }

    @Override
    public void dispose() {
        super.dispose();
        removed = true;
        dataManager = null;
        dataSource = null;
        alignmentTrack = null;
        setVisible(false);
    }

    /**
     * Rescale as necessary, and tell components to repaint
     *
     * @param e
     */
    @Subscribe
    public void receiveDataLoaded(DataLoadedEvent e) {
        ReferenceFrame frame = e.getReferenceFrame();
        //rescale(frame);
        frame.getEventBus().post(new ViewChange.Result());
    }

//    public void rescale(ReferenceFrame iframe) {
//        List<ReferenceFrame> frameList = new ArrayList<ReferenceFrame>();
//        if (iframe != null) frameList.add(iframe);
//        if (globalAutoScale) {
//            frameList.addAll(FrameManager.getFrames());
//        }
//
//        if (autoScale && dataManager != null) {
//
//            int max = 10;
//            for (ReferenceFrame frame : frameList) {
//                AlignmentInterval interval = dataManager.getLoadedInterval(frame.getCurrentRange());
//                if (interval == null) continue;
//
//                int origin = (int) frame.getOrigin();
//                int end = (int) frame.getEnd() + 1;
//
//                int intervalMax = interval.getMaxCount(origin, end);
//                max = intervalMax > max ? intervalMax : max;
//            }
//
//            DataRange newRange = new DataRange(0, max);
//            newRange.setType(getDataRange().getType());
//            super.setDataRange(newRange);
//
//        }
//    }


    public void render(RenderContext context, Rectangle rect) {

        if (context.getScale() > dataManager.getMinVisibleScale() && dataSource == null) {
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

        float maxRange = PreferenceManager.getInstance().getAsFloat(PreferenceManager.SAM_MAX_VISIBLE_RANGE);
        float minVisibleScale = (maxRange * 1000) / 700;

        if (context.getScale() < minVisibleScale && !context.getChr().equals(Globals.CHR_ALL)) {
            //Show coverage calculated from intervals if zoomed in enough
            AlignmentInterval interval = null;
            if (dataManager != null) {
                dataManager.load(context.getReferenceFrame(), renderOptions, true);
                interval = dataManager.getLoadedInterval(context.getReferenceFrame().getCurrentRange());
            }
            if (interval != null) {
                if (interval.contains(context.getChr(), (int) context.getOrigin(), (int) context.getEndLocation())) {
                    //if (autoScale) rescale(context.getReferenceFrame());
                    intervalRenderer.paint(context, rect, interval.getCounts());
                    return;
                }
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
            int startIdx = FeatureUtils.getIndexBefore(start, inViewScores);
            int endIdx = inViewScores.size() - 1;   // Starting guess
            int tmp = FeatureUtils.getIndexBefore(end, inViewScores);
            for (int i = tmp; i < inViewScores.size(); i++) {
                if (inViewScores.get(i).getStart() > end) {
                    endIdx = i - 1;
                    break;
                }
            }
            endIdx = Math.max(startIdx + 1, endIdx);

            return startIdx == 0 && endIdx == inViewScores.size() - 1 ?
                    inViewScores :
                    inViewScores.subList(startIdx, endIdx);
        }
        return inViewScores;
    }


    @Override
    public Range getInViewRange(ReferenceFrame frame) {

        if (dataManager == null || frame.getScale() > dataManager.getMinVisibleScale()) {

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

            AlignmentInterval interval = dataManager.getLoadedInterval(frame.getCurrentRange());
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

    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {

        float maxRange = PreferenceManager.getInstance().getAsFloat(PreferenceManager.SAM_MAX_VISIBLE_RANGE);
        float minVisibleScale = (maxRange * 1000) / 700;
        if (frame.getScale() < minVisibleScale) {
            AlignmentInterval interval = dataManager.getLoadedInterval(frame.getCurrentRange());
            if (interval == null) return null;

            if (interval.contains(chr, (int) position, (int) position)) {
                AlignmentCounts counts = interval.getCounts();
                if (counts != null) {
                    return counts.getValueStringAt((int) position);
                }
            }
        } else {
            return getPrecomputedValueString(chr, position, frame);
        }
        return null;
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

    /**
     * Class to render coverage track, including mismatches.
     * <p/>
     * NOTE:  This class has been extensively optimized with the aid of a profiler,  attempts to "clean up" this code
     * should be done with frequent profiling, or it will likely have detrimental performance impacts.
     */
    class IntervalRenderer {

        private void paint(RenderContext context, Rectangle rect, AlignmentCounts alignmentCounts) {

            Color color = getColor();
            Graphics2D graphics = context.getGraphic2DForColor(color);

            DataRange range = getDataRange();
            double maxRange = range.isLog() ? Math.log10(range.getMaximum() + 1) : range.getMaximum();

            final double rectX = rect.getX();
            final double rectMaxX = rect.getMaxX();
            final double rectY = rect.getY();
            final double rectMaxY = rect.getMaxY();
            final double rectHeight = rect.getHeight();
            final double origin = context.getOrigin();
            final double colorScaleMax = getColorScale().getMaximum();
            final double scale = context.getScale();

            boolean bisulfiteMode = dataManager.getExperimentType() == AlignmentTrack.ExperimentType.BISULFITE;


            // First pass, coverage
            int lastpX = -1;

            int start = alignmentCounts.getStart();
            int step = alignmentCounts.getBucketSize();
            int nPoints = alignmentCounts.getNumberOfPoints();
            boolean isSparse = alignmentCounts instanceof SparseAlignmentCounts;

            for (int idx = 0; idx < nPoints; idx++) {

                int pos = isSparse ? ((SparseAlignmentCounts) alignmentCounts).getPosition(idx) : start + idx * step;
                int pX = (int) (rectX + (pos - origin) / scale);

                if (pX > rectMaxX) {
                    break; // We're done,  data is position sorted so we're beyond the right-side of the view
                }

                int dX = (int) (rectX + (pos + step - origin) / scale) - pX;
                dX = dX < 1 ? 1 : dX;
                if (pX + dX > lastpX) {
                    int pY = (int) rectMaxY - 1;
                    int totalCount = alignmentCounts.getTotalCount(pos);
                    double tmp = range.isLog() ? Math.log10(totalCount + 1) / maxRange : totalCount / maxRange;
                    int height = (int) (tmp * rectHeight);

                    height = Math.min(height, rect.height - 1);
                    int topY = (pY - height);
                    if (dX > 3) {
                        dX--; // Create a little space between bars when there is room.
                    }
                    if (height > 0) {
                        graphics.fillRect(pX, topY, dX, height);

                    }
                    lastpX = pX + dX;
                }
            }

            // Second pass - mark mismatches

            if (alignmentCounts.hasBaseCounts()) {

                lastpX = -1;
                BisulfiteCounts bisulfiteCounts = alignmentCounts.getBisulfiteCounts();
                final int intervalEnd = alignmentCounts.getEnd();
                final int intervalStart = alignmentCounts.getStart();
                byte[] refBases = null;

                // Dont try to compute mismatches for intervals > 10 MB
                if ((intervalEnd - intervalStart) < TEN_MB) {
                    refBases = genome.getSequence(context.getChr(), intervalStart, intervalEnd);
                }


                start = alignmentCounts.getStart();
                nPoints = alignmentCounts.getNumberOfPoints();
                isSparse = alignmentCounts instanceof SparseAlignmentCounts;

                for (int idx = 0; idx < nPoints; idx++) {
                    int pos = isSparse ? ((SparseAlignmentCounts) alignmentCounts).getPosition(idx) : start + idx;

                    BisulfiteCounts.Count bc = null;
                    if (bisulfiteMode && bisulfiteCounts != null) {
                        bc = bisulfiteCounts.getCount(pos);
                    }

                    int pX = (int) (rectX + (pos - origin) / scale);

                    if (pX > rectMaxX) {
                        break; // We're done,  data is position sorted so we're beyond the right-side of the view
                    }

                    int dX = (int) (rectX + (pos + 1 - origin) / scale) - pX;
                    dX = dX < 1 ? 1 : dX;
                    if (pX + dX > lastpX) {


                        // Test to see if any single nucleotide mismatch  (nucleotide other than the reference)
                        // has a quality weight > 20% of the total
                        // Skip this test if the position is in the list of known snps or if the reference is unknown
                        boolean mismatch = false;

                        if (refBases != null) {
                            int refIdx = pos - intervalStart;
                            if (refIdx >= 0 && refIdx < refBases.length) {
                                if (bisulfiteMode) {
                                    mismatch = (bc != null && (bc.methylatedCount + bc.unmethylatedCount) > 0);
                                } else {
                                    byte ref = refBases[refIdx];
                                    mismatch = alignmentCounts.isMismatch(pos, ref, context.getChr(), snpThreshold);
                                }
                            }
                        }

                        if (!mismatch) {
                            continue;
                        }

                        int pY = (int) rectMaxY - 1;

                        int totalCount = alignmentCounts.getTotalCount(pos);
                        double tmp = range.isLog() ? Math.log10(totalCount + 1) / maxRange : totalCount / maxRange;
                        int height = (int) (tmp * rectHeight);

                        height = Math.min(height, rect.height - 1);

                        if (dX > 3) {
                            dX--; // Create a little space between bars when there is room.
                        }

                        if (height > 0) {
                            if (bisulfiteMode) {
                                if (bc != null) {
                                    drawBarBisulfite(context, pos, rect, totalCount, maxRange,
                                            pY, pX, dX, bc, range.isLog());
                                }
                            } else {
                                drawBar(context, pos, rect, totalCount, maxRange,
                                        pY, pX, dX, alignmentCounts, range.isLog());
                            }
                        }
                        lastpX = pX + dX;

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
     * @param rect
     * @param max
     * @param pY
     * @param pX
     * @param dX
     * @param interval
     * @return
     */

    int drawBar(RenderContext context,
                int pos,
                Rectangle rect,
                double totalCount,
                double max,
                int pY,
                int pX,
                int dX,
                AlignmentCounts interval,
                boolean isLog) {

        for (char nucleotide : nucleotides) {

            int count = interval.getCount(pos, (byte) nucleotide);

            Color c = AlignmentRenderer.nucleotideColors.get(nucleotide);

            Graphics2D tGraphics = context.getGraphic2DForColor(c);

            double tmp = isLog ?
                    (count / totalCount) * Math.log10(totalCount + 1) / max :
                    count / max;
            int height = (int) (tmp * rect.getHeight());

            height = Math.min(pY - rect.y, height);
            int baseY = pY - height;

            if (height > 0) {
                tGraphics.fillRect(pX, baseY, dX, height);
            }

            pY = baseY;
        }
        return pX + dX;
    }

    int drawBarBisulfite(RenderContext context,
                         int pos,
                         Rectangle rect,
                         double totalCount,
                         double maxRange,
                         int pY,
                         int pX0,
                         int dX,
                         BisulfiteCounts.Count count,
                         boolean isLog) {

        // If bisulfite mode, we expand the rectangle to make it more visible.  This code is copied from
        // AlignmentRenderer
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

        double tmp = isLog ?
                (nMethylated / totalCount) * Math.log10(totalCount + 1) / maxRange :
                nMethylated / maxRange;
        int height = (int) (tmp * rect.getHeight());

        height = Math.min(pY - rect.y, height);
        int baseY = pY - height;
        if (height > 0) {
            tGraphics.fillRect(pX, baseY, dX, height);
        }
        pY = baseY;

        c = Color.blue;
        tGraphics = context.getGraphic2DForColor(c);

        tmp = isLog ?
                (unMethylated / totalCount) * Math.log10(totalCount + 1) / maxRange :
                unMethylated / maxRange;
        height = (int) (tmp * rect.getHeight());

        height = Math.min(pY - rect.y, height);
        baseY = pY - height;
        if (height > 0) {
            tGraphics.fillRect(pX, baseY, dX, height);

        }
        return pX + dX;
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

            Color c = AlignmentRenderer.nucleotideColors.get(nucleotide);
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
        int minQ = prefs.getAsInt(PreferenceManager.SAM_BASE_QUALITY_MIN);
        ColorUtilities.getRGBColorComponents(color);
        if (qual < minQ) {
            alpha = 0.1f;
        } else {
            int maxQ = prefs.getAsInt(PreferenceManager.SAM_BASE_QUALITY_MAX);
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

        popupMenu.addSeparator();
        addShowItems(popupMenu);

        return popupMenu;
    }

    private void addCopyDetailsItem(IGVPopupMenu popupMenu, TrackClickEvent te) {
        JMenuItem copyDetails = new JMenuItem("Copy Details to Clipboard");
        copyDetails.setEnabled(false);
        if (te.getFrame() != null) {
            final String details = getValueStringAt(te.getFrame().getChrName(), te.getChromosomePosition(), te.getMouseEvent().getY(), te.getFrame());
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

        maxValItem.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {

                String value = JOptionPane.showInputDialog("Allele frequency threshold: ", Float.valueOf(snpThreshold));
                if (value == null) {
                    return;
                }
                try {
                    float tmp = Float.parseFloat(value);
                    snpThreshold = tmp;
                    IGV.getInstance().repaintDataPanels();
                } catch (Exception exc) {
                    //log
                }

            }
        });
        menu.add(maxValItem);

        return maxValItem;
    }

    public void addShowItems(JPopupMenu menu) {

        if (alignmentTrack != null) {
            final JMenuItem alignmentItem = new JCheckBoxMenuItem("Show Alignment Track");
            alignmentItem.setSelected(alignmentTrack.isVisible());
            alignmentItem.setEnabled(!alignmentTrack.isRemoved());
            alignmentItem.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    alignmentTrack.onAlignmentTrackEvent(new AlignmentTrackEvent(CoverageTrack.this, AlignmentTrackEvent.Type.VISIBLE, alignmentItem.isSelected()));
                    if (alignmentItem.isSelected()) {
                        alignmentTrack.onAlignmentTrackEvent(new AlignmentTrackEvent(CoverageTrack.this, AlignmentTrackEvent.Type.RELOAD));
                    }
                }
            });
            menu.add(alignmentItem);

            final SpliceJunctionTrack spliceJunctionTrack = alignmentTrack.getSpliceJunctionTrack();
            if (spliceJunctionTrack != null) {
                final JMenuItem junctionItem = new JCheckBoxMenuItem("Show Splice Junction Track");
                junctionItem.setSelected(spliceJunctionTrack.isVisible());
                junctionItem.setEnabled(!spliceJunctionTrack.isRemoved());
                junctionItem.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        alignmentTrack.onAlignmentTrackEvent(new AlignmentTrackEvent(CoverageTrack.this, AlignmentTrackEvent.Type.SPLICE_JUNCTION, junctionItem.isSelected()));
                        IGV.getInstance().getMainPanel().revalidate();
                    }
                });

                menu.add(junctionItem);
            }

            final JMenuItem coverageItem = new JMenuItem("Hide Track");
            coverageItem.setEnabled(!isRemoved());
            coverageItem.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    UIUtilities.invokeOnEventThread(new Runnable() {

                        public void run() {
                            setVisible(false);
                            IGV.getInstance().getMainPanel().revalidate();

                        }
                    });
                }
            });
            menu.add(coverageItem);
        }
    }


    public void addLoadCoverageDataItem(JPopupMenu menu) {
        // Change track height by attribute
        final JMenuItem item = new JCheckBoxMenuItem("Load pre-computed coverage data...");
        item.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {

                final PreferenceManager prefs = PreferenceManager.getInstance();
                File initDirectory = prefs.getLastTrackDirectory();
                File file = FileDialogUtils.chooseFile("Select coverage file", initDirectory, FileDialog.LOAD);
                if (file != null) {
                    prefs.setLastTrackDirectory(file.getParentFile());
                    String path = file.getAbsolutePath();
                    if (path.endsWith(".tdf") || path.endsWith(".tdf")) {
                        TDFReader reader = TDFReader.getReader(file.getAbsolutePath());
                        TDFDataSource ds = new TDFDataSource(reader, 0, getName() + " coverage", genome);
                        setDataSource(ds);
                        IGV.getInstance().repaintDataPanels();
                    } else if (path.endsWith(".counts")) {
                        CoverageDataSource ds = new GobyCountArchiveDataSource(file);
                        setDataSource(ds);
                        IGV.getInstance().repaintDataPanels();
                    } else {
                        MessageUtils.showMessage("Coverage data must be in .tdf format");
                    }
                }
            }
        });

        item.setEnabled(dataSource == null);
        menu.add(item);

    }


    @SubtlyImportant
    private static CoverageTrack getNextTrack() {
        return (CoverageTrack) IGVSessionReader.getNextTrack();
    }

    public void setGlobalAutoScale(boolean globalAutoScale) {
        this.globalAutoScale = globalAutoScale;
    }

}
