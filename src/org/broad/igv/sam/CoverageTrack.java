/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam;

import com.jidesoft.swing.JidePopupMenu;
import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ColorUtilities;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.renderer.BarChartRenderer;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.DataRenderer;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.tdf.TDFDataSource;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.track.*;
import org.broad.igv.ui.DataRangeDialog;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGVMainFrame;
import org.broad.igv.ui.util.FileChooserDialog;
import org.broad.igv.ui.util.MessageUtils;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class CoverageTrack extends AbstractTrack {

    private static Logger log = Logger.getLogger(CoverageTrack.class);
    private static Map<String, Set<Integer>> knownSnps;

    char[] nucleotides = {'a', 'c', 'g', 't', 'n'};
    public static Color lightBlue = new Color(0, 0, 150);
    private float[] bgColorComps = new float[3];
    private static final boolean DEFAULT_AUTOSCALE = true;
    private static final boolean DEFAULT_SHOW_REFERENCE = false;

    // User settable state -- these attributes should be stored in the session file
    boolean showReference;
    boolean autoScale = DEFAULT_AUTOSCALE;
    private float snpThreshold;
    private boolean showAllSnps;

    AlignmentDataManager dataManager;
    TDFDataSource dataSource;
    DataRenderer dataSourceRenderer; // = new BarChartRenderer();
    IntervalRenderer intervalRenderer;
    PreferenceManager prefs;
    JMenuItem dataRangeItem;
    JMenuItem autoscaleItem;


    public CoverageTrack(String id, String name) {
        super(id, name);
        super.setDataRange(new DataRange(0, 0, 60));
        intervalRenderer = new IntervalRenderer();

        setColor(AlignmentRenderer.grey1);
        AlignmentRenderer.grey1.getColorComponents(bgColorComps);

        prefs = PreferenceManager.getInstance();
        snpThreshold = prefs.getAsFloat(PreferenceManager.SAM_ALLELE_THRESHOLD);
        autoScale = DEFAULT_AUTOSCALE;
        showReference = DEFAULT_SHOW_REFERENCE;
        showAllSnps = prefs.getAsBoolean(PreferenceManager.COVERAGE_SHOW_ALL_MISMATCHES);
        //TODO  logScale = prefs.

        String snpsFile = prefs.get(PreferenceManager.KNOWN_SNPS, null);
        if (snpsFile != null && knownSnps == null) {
            loadKnownSnps(snpsFile);
        }
    }

    public void setDataManager(AlignmentDataManager dataManager) {
        this.dataManager = dataManager;
    }

    public void setDataSource(TDFDataSource dataSource) {
        this.dataSource = dataSource;
        dataSourceRenderer = new BarChartRenderer();
        setDataRange(new DataRange(0, 0, 1.5f * (float) dataSource.getDataMax()));

    }


    @Override
    public void setDataRange(DataRange axisDefinition) {
        // Explicitly setting a data range turns off auto-scale
        autoScale = false;
        super.setDataRange(axisDefinition);
    }

    public void rescale() {
        if (autoScale & dataManager != null) {
            for (AlignmentInterval interval : dataManager.getLoadedIntervals()) {
                rescaleInterval(interval);
            }
        }
    }

    public void rescale(ReferenceFrame frame) {
        rescaleInterval(dataManager.getLoadedInterval(frame));
    }

    private void rescaleInterval(AlignmentInterval interval) {
        if (interval != null) {
            int max = Math.max(10, interval.getMaxCount());
            DataRange.Type type = getDataRange().getType();
            super.setDataRange(new DataRange(0, 0, max));
            getDataRange().setType(type);
        }
    }


    public void render(RenderContext context, Rectangle rect) {

        float maxRange = PreferenceManager.getInstance().getAsFloat(PreferenceManager.SAM_MAX_VISIBLE_RANGE);
        float minVisibleScale = (maxRange * 1000) / 700;

        if (context.getScale() < minVisibleScale) {

            AlignmentInterval interval = null;
            if (dataManager != null) {
                interval = dataManager.getLoadedInterval(context);
                
            }




            if (interval != null && interval.contains(context.getGenomeId(), context.getChr(), (int) context.getOrigin(),
                    (int) context.getEndLocation())) {
                List<AlignmentCounts> counts = interval.counts;
                intervalRenderer.paint(context, rect, counts);
            }
        }
        // Use precomputed data source, if any
        else if (dataSource != null) {
            String chr = context.getChr();
            int start = (int) context.getOrigin();
            int end = (int) context.getEndLocation();
            int zoom = context.getZoom();
            List<LocusScore> scores = dataSource.getSummaryScoresForRange(chr, start, end, zoom);
            if (scores != null) {
                dataSourceRenderer.render(scores, context, rect, this);
            }

        }
        drawBorder(context, rect);


    }

    private void drawBorder(RenderContext context, Rectangle rect) {
        // Draw border
        context.getGraphic2DForColor(Color.gray).drawLine(
                rect.x, rect.y + rect.height,
                rect.x + rect.width, rect.y + rect.height);

        // Draw scale
        DataRange range = getDataRange();
        if (range != null) {
            Graphics2D g = context.getGraphic2DForColor(Color.black);
            Font font = g.getFont();
            Font smallFont = FontManager.getScalableFont(8);
            try {
                g.setFont(smallFont);
                String scale = "[" + (int) range.getMinimum() + " - " +
                        (int) range.getMaximum() + "]";
                g.drawString(scale, rect.x + 5, rect.y + 10);

            } finally {
                g.setFont(font);
            }
        }
    }

    public void setStatType(WindowFunction type) {
    }

    public WindowFunction getWindowFunction() {
        return null;
    }

    public void setRendererClass(Class rc) {
    }

    public Renderer getRenderer() {
        return null;
    }

    public boolean isLogNormalized() {
        return false;
    }

    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {

        float maxRange = PreferenceManager.getInstance().getAsFloat(PreferenceManager.SAM_MAX_VISIBLE_RANGE);
        float minVisibleScale = (maxRange * 1000) / 700;
        if (frame.getScale() < minVisibleScale) {
            AlignmentInterval interval = dataManager.getLoadedInterval(frame);
            if (interval != null && interval.contains(chr, (int) position, (int) position)) {
                StringBuffer buf = new StringBuffer();
                final int pos = (int) position - 1;
                int totalCount = interval.getTotalCount(pos);
                buf.append("Total count: " + totalCount + "<br>");
                for (char c : nucleotides) {
                    int negCount = interval.getNegCount(pos, (byte) c);
                    int posCount = interval.getPosCount(pos, (byte) c);
                    int count = negCount + posCount;
                    int percent = (int) Math.round(((float) count) * 100 / totalCount);
                    char cU = Character.toUpperCase(c);
                    buf.append(cU + "      : " + count);
                    if (count == 0) {
                        buf.append("<br>");
                    } else {
                        buf.append("  (" + percent + "%,     " + posCount + "+,   " + negCount + "- )<br>");
                    }
                }
                return buf.toString();
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
            LocusScore score = (LocusScore) FeatureUtils.getFeatureAt(position, minWidth, scores);
            return score == null ? "" : "Mean count: " + score.getScore();
        }
    }

    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, ReferenceFrame frame) {
        return 0;
    }


    /**
     * Load the set of known snps from a tab delimited file, format
     * chr < tab> location
     * The location is "1 base"  (first nucleotide is position 1).
     *
     * @param snpFile
     */
    private static synchronized void loadKnownSnps(String snpFile) {

        // This method might get called many times concurrently, but we only want to load these once.
        if (knownSnps != null) {
            return;
        }

        knownSnps = new HashMap();
        AsciiLineReader reader = null;
        try {
            reader = ParsingUtils.openAsciiReader(new ResourceLocator(snpFile));
            String nextLine = "";
            while ((nextLine = reader.readLine()) != null) {
                String[] tokens = nextLine.split("\t");
                String chr = tokens[0];
                Set<Integer> snps = knownSnps.get(chr);
                if (snps == null) {
                    snps = new HashSet(10000);
                    knownSnps.put(chr, snps);
                }
                snps.add(new Integer(tokens[1]));
            }
        } catch (Exception e) {
            knownSnps = null;
            log.error("", e);
            MessageUtils.showMessage("Error loading snps file: " + snpFile + " (" + e.toString() + ")");
        } finally {
            reader.close();
        }


    }

    /*
                           if (plusMinusOption || nonRefOption) {

                            int pY = (int) (rect.getY() + rect.getMaxY()) / 2;

                            int totalNegCount = interval.getNegTotal(pos);
                            int height = (int) (totalNegCount * rect.getHeight() / getColorScale().getMaximum());
                            height = Math.min(height, rect.height / 2 - 1);
                            if (height > 0) {

                                if (!nonRefOption) {
                                    graphics.fillRect(pX, pY, dX, height);

                                    if (colorBases) {
                                        for (char c : nucleotides) {
                                            if (nonRefOption == false || c != ref) {
                                                pY = drawBar(context, pos, rect, getColorScale().getMaximum(), pY, pX, dX, c,
                                                        plusMinusOption || nonRefOption, false, interval);
                                            }
                                        }
                                    }
                                }
                            }


                            pY = (int) (rect.getY() + rect.getMaxY()) / 2;
                            int totalPosCount = interval.getPosTotal(pos);

                            height = (int) (totalPosCount * rect.getHeight() / getColorScale().getMaximum());
                            height = Math.min(height, rect.height / 2 - 1);
                            int topY = (pY - height);

                            if (height > 0) {

                                if (!nonRefOption) {
                                    graphics.fillRect(pX, topY, dX, height);
                                }

                                if (colorBases) {
                                    for (char c : nucleotides) {
                                        if (nonRefOption == false || c != ref) {
                                            pY = drawBar(context, pos, rect, getColorScale().getMaximum(), pY, pX, dX, c,
                                                    plusMinusOption || nonRefOption, true, interval);
                                        }
                                    }
                                }
                            }

                            pY = (int) (rect.getY() + rect.getMaxY()) / 2;
                            Graphics2D blackGraphics = context.getGraphic2DForColor(lightBlue);
                            blackGraphics.drawLine(0, pY, rect.width, pY);


                        } else {

                            int totalCount = interval.getTotalCount(pos);

                            int pY = (int) rect.getMaxY() - 1;
                            int height = (int) (totalCount * rect.getHeight() / getColorScale().getMaximum());
                            height = Math.min(height, rect.height - 1);
                            int topY = (pY - height);

                            if (height > 0) {
                                graphics.fillRect(pX, topY, dX, height);

                                if (colorBases) {
                                    for (char c : nucleotides) {
                                        pY = drawBar(context, pos, rect, dataMax, pY, pX, dX, c, plusMinusOption, true, interval);
                                    }
                                }
                            }
                        }
     */

    class IntervalRenderer {


        private void paint(RenderContext context, Rectangle rect, List<AlignmentCounts> countList) {

            Graphics2D graphics = context.getGraphic2DForColor(AlignmentRenderer.grey1);

            DataRange range = getDataRange();
            double max = range.isLog() ? Math.log10(range.getMaximum()) : range.getMaximum();

            // Temporary until proper windowing is implemented
            int lastpX = -1;

            for (AlignmentCounts alignmentCounts : countList) {

                for (int pos = alignmentCounts.getStart(); pos < alignmentCounts.getEnd(); pos++) {

                    int pX = (int) (rect.getX() + (pos - context.getOrigin()) / context.getScale());
                    int dX = Math.max(1,
                            (int) (rect.getX() + (pos + 1 - context.getOrigin()) / context.getScale()) - pX);
                    if (dX > 3) {
                        dX--;
                    }

                    if (pX > rect.getMaxX()) {
                        break;
                    }

                    if (pX + dX >= 0) {

                        // Test to see if any single nucleotide mismatch  (nucleotide other than the reference)
                        // has posA quality weight > 20% of the total
                        // Skip this test if the position is in the list of known snps
                        char ref = Character.toLowerCase((char) alignmentCounts.getReference(pos));
                        boolean mismatch = false;

                        Set<Integer> snps = knownSnps == null ? null : knownSnps.get(context.getChr());
                        if (snps == null || !snps.contains(pos + 1)) {
                            float threshold = snpThreshold * alignmentCounts.getTotalQuality(pos);
                            if (ref > 0) {
                                for (char c : nucleotides) {
                                    if (c != ref && c != 'n' && alignmentCounts.getQuality(pos, (byte) c) > threshold) {
                                        mismatch = true;
                                        break;
                                    }
                                }
                            }
                        }

                        if (pX > lastpX || mismatch) {

                            boolean strandOption = false;

                            if (strandOption) {

                                int pY = (int) (rect.getY() + rect.getMaxY()) / 2;

                                int totalNegCount = alignmentCounts.getNegTotal(pos);
                                int height = (int) (totalNegCount * rect.getHeight() / getColorScale().getMaximum());
                                height = Math.min(height, rect.height / 2 - 1);
                                if (height > 0) {

                                    graphics.fillRect(pX, pY, dX, height);

                                    if (mismatch || showAllSnps) {
                                        for (char c : nucleotides) {
                                            if (c != ref) {
                                                pY = drawStrandBar(context, pos, rect, getColorScale().getMaximum(), pY, pX, dX, c,
                                                        false, alignmentCounts);
                                            }
                                        }
                                    }

                                }

                                pY = (int) (rect.getY() + rect.getMaxY()) / 2;
                                int totalPosCount = alignmentCounts.getPosTotal(pos);

                                height = (int) (totalPosCount * rect.getHeight() / getColorScale().getMaximum());
                                height = Math.min(height, rect.height / 2 - 1);
                                int topY = (pY - height);

                                if (height > 0) {

                                    graphics.fillRect(pX, topY, dX, height);

                                    if (mismatch || showAllSnps) {
                                        for (char c : nucleotides) {
                                            if (c != ref) {
                                                pY = drawStrandBar(context, pos, rect, getColorScale().getMaximum(), pY, pX, dX, c,
                                                        true, alignmentCounts);
                                            }
                                        }
                                    }
                                }

                                pY = (int) (rect.getY() + rect.getMaxY()) / 2;
                                Graphics2D blackGraphics = context.getGraphic2DForColor(lightBlue);
                                blackGraphics.drawLine(0, pY, rect.width, pY);
                            } else {

                                int pY = (int) rect.getMaxY() - 1;

                                int totalCount = alignmentCounts.getTotalCount(pos);

                                double tmp = range.isLog() ? Math.log10(totalCount) / max : totalCount / max;
                                int height = (int) (tmp * rect.getHeight());

                                height = Math.min(height, rect.height - 1);
                                int topY = (pY - height);

                                if (height > 0) {
                                    graphics.fillRect(pX, topY, dX, height);

                                    if (mismatch || showAllSnps) {
                                        for (char c : nucleotides) {
                                            if (!(c == ref && showAllSnps)) {
                                                pY = drawBar(context, pos, rect, totalCount, max,
                                                        pY, pX, dX, c, alignmentCounts, range.isLog());
                                            }
                                        }
                                    }
                                }
                            }
                            lastpX = pX;
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
         * @param nucleotide
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
                    char nucleotide,
                    AlignmentCounts interval,
                    boolean isLog) {

            int count = interval.getCount(pos, (byte) nucleotide);

            Color c = AlignmentRenderer.getNucleotideColors().get(nucleotide);

            if (showAllSnps) {
                int q = interval.getAvgQuality(pos, (byte) nucleotide);
                c = getShadedColor(q, c);
            }

            Graphics2D tGraphics = context.getGraphic2DForColor(c);

            double tmp = isLog ?
                    (count / totalCount) * Math.log10(totalCount) / max :
                    count / max;
            int height = (int) (tmp * rect.getHeight());

            height = Math.min(pY - rect.y, height);
            int baseY = pY - height;

            if (height > 0) {
                tGraphics.fillRect(pX, baseY, dX, height);
            }
            return baseY;
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
         * @param nucleotide
         * @param isPositive
         * @param interval
         * @return
         */
        int drawStrandBar(RenderContext context,
                          int pos,
                          Rectangle rect,
                          double maxCount,
                          int pY,
                          int pX,
                          int dX,
                          char nucleotide,
                          boolean isPositive,
                          AlignmentCounts interval) {


            Color c = AlignmentRenderer.getNucleotideColors().get(nucleotide);
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
            return isPositive ? baseY : baseY + height;
        }

    }


    static float[] colorComps = new float[3];

    private Color getShadedColor(int qual, Color color) {
        float alpha = 0;
        int minQ = prefs.getAsInt(PreferenceManager.SAM_BASE_QUALITY_MIN);
        color.getRGBColorComponents(colorComps);
        if (qual < minQ) {
            alpha = 0.1f;
        } else {
            int maxQ = prefs.getAsInt(PreferenceManager.SAM_BASE_QUALITY_MAX);
            alpha = Math.max(0.1f, Math.min(1.0f, 0.1f + 0.9f * (qual - minQ) / (maxQ - minQ)));
        }

        // Round alpha to nearest 0.1, for effeciency;
        alpha = ((int) (alpha * 10 + 0.5f)) / 10.0f;
        color = ColorUtilities.getCompositeColor(bgColorComps, colorComps, alpha);
        return color;
    }


    /**
     * Called by session writer.  Return instance variable values as a map of strings.  Used to record current state
     * of object.   Variables with default values are not stored, as it is presumed the user has not changed them.
     *
     * @return
     */
    @Override
    public Map<String, String> getPersistentState() {
        Map<String, String> attributes = super.getPersistentState();
        if (dataSource != null) {
            attributes.put("path", dataSource.getPath());
        }
        prefs = PreferenceManager.getInstance();
        if (snpThreshold != prefs.getAsFloat(PreferenceManager.SAM_ALLELE_THRESHOLD)) {
            attributes.put("snpThreshold", String.valueOf(snpThreshold));
        }
        if (autoScale != DEFAULT_AUTOSCALE) {
            attributes.put("autoScale", String.valueOf(autoScale));
        }
        if (showReference != DEFAULT_SHOW_REFERENCE) {
            attributes.put("showReference", String.valueOf(showReference));
        }
        if (showAllSnps != prefs.getAsBoolean(PreferenceManager.COVERAGE_SHOW_ALL_MISMATCHES)) {
            attributes.put("showAllSnps", String.valueOf(showAllSnps));
        }

        return attributes;
    }

    /**
     * Called by session reader.  Restores state of object.
     *
     * @param attributes
     */
    @Override
    public void restorePersistentState(Map<String, String> attributes) {
        super.restorePersistentState(attributes);    //To change body of overridden methods use File | Settings | File Templates.

        String value;
        value = attributes.get("path");
        if (value != null) {
            TDFReader reader = TDFReader.getReader(value);
            TDFDataSource ds = new TDFDataSource(reader, 0, "");
            setDataSource(ds);
        }
        value = attributes.get("snpThreshold");
        if (value != null) {
            snpThreshold = Float.parseFloat(value);
        }
        value = attributes.get("autoScale");
        if (value != null) {
            autoScale = Boolean.parseBoolean(value);
        }
        value = attributes.get("showReference");
        if (value != null) {
            showReference = Boolean.parseBoolean(value);
        }
        value = attributes.get("showAllSnps");
        if (value != null) {
            showAllSnps = Boolean.parseBoolean(value);
        }

    }

    /**
     * Override to return a specialized popup menu
     *
     * @return
     */
    @Override
    public JPopupMenu getPopupMenu(TrackClickEvent te) {

        JPopupMenu popupMenu = new JidePopupMenu();

        JLabel popupTitle = new JLabel("  " + getName(), JLabel.CENTER);

        Font newFont = popupMenu.getFont().deriveFont(Font.BOLD, 12);
        popupTitle.setFont(newFont);
        if (popupTitle != null) {
            popupMenu.add(popupTitle);
        }

        popupMenu.addSeparator();

        // addSortMenuItem(popupMenu);
        // addPackMenuItem(popupMenu);
        // addShadeBaseMenuItem(popupMenu);
        // addCopyToClipboardItem(popupMenu, evt);
        // addGoToMate(popupMenu, evt);
        // popupMenu.addSeparator();


        //JLabel trackSettingsHeading = new JLabel("  Track Settings",
        //        JLabel.LEFT);
        //trackSettingsHeading.setFont(newFont);

        //popupMenu.add(trackSettingsHeading);

        ArrayList<Track> tmp = new ArrayList();
        tmp.add(this);
        popupMenu.add(TrackMenuUtils.getTrackRenameItem(tmp));


        addAutoscaleItem(popupMenu);
        addLogScaleItem(popupMenu);
        dataRangeItem = addDataRangeItem(popupMenu, tmp);
        dataRangeItem.setEnabled(!autoScale);

        this.addSnpTresholdItem(popupMenu);

        popupMenu.addSeparator();
        addLoadCoverageDataItem(popupMenu);
        popupMenu.addSeparator();

        popupMenu.add(TrackMenuUtils.getRemoveMenuItem(tmp));

        return popupMenu;
    }


    public JMenuItem addDataRangeItem(JPopupMenu menu, final Collection<Track> selectedTracks) {
        JMenuItem maxValItem = new JMenuItem("Set Data Range");

        maxValItem.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                if (selectedTracks.size() > 0) {

                    DataRange prevAxisDefinition = selectedTracks.iterator().next().getDataRange();
                    DataRangeDialog dlg = new DataRangeDialog(
                            IGVMainFrame.getInstance(),
                            prevAxisDefinition);
                    dlg.setHideMid(true);
                    dlg.setVisible(true);
                    if (!dlg.isCanceled()) {
                        float min = Math.min(dlg.getMin(), dlg.getMax());
                        float max = Math.max(dlg.getMin(), dlg.getMax());
                        float mid = dlg.getBase();
                        if (mid < min) mid = min;
                        else if (mid > max) mid = max;
                        DataRange dataRange = new DataRange(min, mid, max);
                        dataRange.setType(getDataRange().getType());

                        // dlg.isFlipAxis());
                        for (Track track : selectedTracks) {
                            track.setDataRange(dataRange);
                        }
                        IGVMainFrame.getInstance().repaint();
                    }
                }

            }
        });
        menu.add(maxValItem);

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
                    IGVMainFrame.getInstance().repaintDataPanels();
                }
                catch (Exception exc) {
                    //log
                }

            }
        });
        menu.add(maxValItem);

        return maxValItem;
    }

    public void addLoadCoverageDataItem(JPopupMenu menu) {
        // Change track height by attribute
        final JMenuItem item = new JCheckBoxMenuItem("Load coverage data...");
        item.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {


                FileChooserDialog trackFileDialog = IGVMainFrame.getInstance().getTrackFileChooser();
                trackFileDialog.setMultiSelectionEnabled(false);
                trackFileDialog.setVisible(true);
                if (!trackFileDialog.isCanceled()) {
                    File file = trackFileDialog.getSelectedFile();
                    String path = file.getAbsolutePath();
                    if (path.endsWith(".tdf") || path.endsWith(".tdf")) {

                        TDFReader reader = TDFReader.getReader(file.getAbsolutePath());
                        TDFDataSource ds = new TDFDataSource(reader, 0, getName() + " coverage");
                        setDataSource(ds);
                        IGVMainFrame.getInstance().repaintDataPanels();
                    } else {
                        MessageUtils.showMessage("Coverage data must be in .tdf format");
                    }
                }
            }
        });

        menu.add(item);
    }

    public void addAutoscaleItem(JPopupMenu menu) {
        // Change track height by attribute
        autoscaleItem = new JCheckBoxMenuItem("Autoscale");
        autoscaleItem.setSelected(autoScale);
        autoscaleItem.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {

                autoScale = autoscaleItem.isSelected();
                dataRangeItem.setEnabled(!autoScale);
                if (autoScale) {
                    rescale();
                }
                IGVMainFrame.getInstance().repaintDataPanels();

            }
        });

        menu.add(autoscaleItem);
    }

    public void addLogScaleItem(JPopupMenu menu) {
        // Change track height by attribute
        final DataRange dataRange = getDataRange();
        final JCheckBoxMenuItem logScaleItem = new JCheckBoxMenuItem("Log scale");
        final boolean logScale = dataRange.getType() == DataRange.Type.LOG;
        logScaleItem.setSelected(logScale);
        logScaleItem.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {

                DataRange.Type scaleType = logScaleItem.isSelected() ?
                        DataRange.Type.LOG :
                        DataRange.Type.LINEAR;
                dataRange.setType(scaleType);
                IGVMainFrame.getInstance().repaintDataPanels();
            }
        });

        menu.add(logScaleItem);
    }

}
