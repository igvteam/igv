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

package org.broad.igv.peaks;

import org.broad.igv.data.LocusScoreUtils;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.tdf.TDFDataSource;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.*;
import java.io.IOException;
import java.lang.ref.SoftReference;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 * @date Apr 22, 2011
 */
public class PeakTrack extends AbstractTrack {


    static List<SoftReference<PeakTrack>> instances = new ArrayList();

    private static PeakControlDialog controlDialog;
    private static float scoreThreshold = 30;
    private static float foldChangeThreshold = 0;

    private static ColorOption colorOption = ColorOption.SCORE;
    private static boolean showPeaks = true;
    private static boolean showSignals = false;

    int nTimePoints;
    Map<String, List<Peak>> peakMap = new HashMap();
    Map<String, List<Peak>> filteredPeakMap = new HashMap();
    Renderer renderer = new PeakRenderer();

    // Path the signal (TDF) file
    String signalPath;
    TDFDataSource signalSource;

    String[] timeSignalPaths;
    TDFDataSource[] timeSignalSources;

    // Data range
    DataRange scoreDataRange = new DataRange(0, 0, 100);
    DataRange signalDataRange = new DataRange(0, 0, 1000f);

    static boolean commandBarAdded = false;

    int bandHeight;
    int signalHeight;
    int peakHeight;
    int gapHeight;


    public PeakTrack(ResourceLocator locator, Genome genome) throws IOException {
        super(locator);
        setHeight(30);
        loadPeaks(locator.getPath());

        instances.add(new SoftReference(this));

        if (!commandBarAdded) {
            IGV.getInstance().getContentPane().add(new PeakCommandBar());
            commandBarAdded = true;
        }
    }

    private void loadPeaks(String path) throws IOException {

        PeakParser parser = new PeakParser();
        List<Peak> peaks = parser.loadPeaks(path);

        nTimePoints = parser.getnTimePoints();
        signalPath = parser.getSignalPath();
        if (signalPath != null) {
            signalSource = new TDFDataSource(TDFReader.getReader(signalPath), 0, "");
            signalSource.setNormalizeCounts(true, 1.0e9f);
        }

        timeSignalPaths = parser.timeSignalPaths;
        if (timeSignalPaths != null && timeSignalPaths.length > 0) {
            timeSignalSources = new TDFDataSource[timeSignalPaths.length];
            for (int i = 0; i < timeSignalPaths.length; i++) {
                try {
                    timeSignalSources[i] = new TDFDataSource(TDFReader.getReader(timeSignalPaths[i]), 0, "");
                    timeSignalSources[i].setNormalizeCounts(true, 1.0e9f);
                } catch (Exception e) {
                    timeSignalSources[i] = null;
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }
        }

        TrackProperties props = parser.getTrackProperties();
        if (props != null) {
            setProperties(props);
        }

        for (Peak peak : peaks) {
            String chr = peak.getChr();
            List<Peak> peakList = peakMap.get(chr);
            if (peakList == null) {
                peakList = new ArrayList();
                peakMap.put(chr, peakList);
            }
            peakList.add(peak);
        }
    }


    @Override
    public JPopupMenu getPopupMenu(TrackClickEvent te) {
        return new PeakTrackMenu(this);
    }

    @Override
    public DataRange getDataRange() {
        return showSignals ? signalDataRange : scoreDataRange;
    }

    public void render(RenderContext context, Rectangle rect) {

        List<Peak> peakList = getFilteredPeaks(context.getChr());
        if (peakList == null) {
            return;
        }

        renderer.render(peakList, context, rect, this);
    }

    public Renderer getRenderer() {
        return renderer;
    }


    @Override
    public int getMinimumHeight() {
        int h = 0;
        if (showPeaks) h += 5;
        if (showSignals) h += 10;
        if (showPeaks && showSignals) h += 2;

        if (getDisplayMode() == Track.DisplayMode.COLLAPSED) {
            return h;
        } else {
            return (nTimePoints + 1) * h + 2;
        }
    }


    @Override
    public void setHeight(int h) {
        super.setHeight(h);

        int nBands = getDisplayMode() == DisplayMode.COLLAPSED ? 1 : nTimePoints + 1;

        bandHeight = h / nBands;
        peakHeight = Math.max(5, Math.min(bandHeight / 3, 10));
        signalHeight = bandHeight - peakHeight - gapHeight;

    }

    @Override
    public void setDisplayMode(DisplayMode mode) {
        super.setDisplayMode(mode);
        if (mode == Track.DisplayMode.COLLAPSED) {
            setHeight(bandHeight);
        } else {
            setHeight((nTimePoints + 1) * bandHeight + gapHeight);
        }
    }


    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {
        StringBuffer buf = new StringBuffer();
        buf.append(getName());
        if (showPeaks) {
            List<Peak> scores = getFilteredPeaks(chr);
            LocusScore score = getLocusScoreAt(scores, position, frame);
            buf.append((score == null) ? "" : score.getValueString(position, getWindowFunction()));
            if (showSignals) {
                buf.append("<br>");
            }
        }
        if (showSignals && signalSource != null) {
            List<LocusScore> scores = signalSource.getSummaryScoresForRange(chr, (int) frame.getOrigin(), (int) frame.getEnd(), frame.getZoom());
            LocusScore score = getLocusScoreAt(scores, position, frame);
            buf.append((score == null) ? "" : "Score = " + score.getScore());
        }


        return buf.toString();
    }


    // TODO -- the code below is an exact copy of code in DataTrack.   Refactor to share this.

    private LocusScore getLocusScoreAt(List<? extends LocusScore> scores, double position, ReferenceFrame frame) {

        if (scores == null) {
            return null;
        } else {
            // give a 2 pixel window, otherwise very narrow features will be missed.
            double bpPerPixel = frame.getScale();
            int buffer = (int) (2 * bpPerPixel);    /* * */
            return LocusScoreUtils.getFeatureAt(position, buffer, scores);
        }
    }

    public synchronized List<Peak> getFilteredPeaks(String chr) {
        List<Peak> filteredPeaks = filteredPeakMap.get(chr);
        if (filteredPeaks == null) {
            List<Peak> allPeaks = peakMap.get(chr);
            if (allPeaks == null) {
                return null;
            }
            filteredPeaks = new ArrayList(allPeaks.size() / 2);
            for (Peak peak : allPeaks) {
                if (peak.getCombinedScore() >= scoreThreshold &&
                        peak.getFoldChange() >= foldChangeThreshold) {
                    filteredPeaks.add(peak);
                }
            }
            filteredPeakMap.put(chr, filteredPeaks);
        }
        return filteredPeaks;
    }


    private static void clearFilteredLists() {
        for (SoftReference<PeakTrack> instance : instances) {
            PeakTrack track = instance.get();
            if (track != null) {
                track.filteredPeakMap.clear();
            }
        }
    }


    public static boolean controlDialogIsOpen() {
        return controlDialog != null && controlDialog.isVisible();
    }


    static synchronized void openControlDialog() {
        if (controlDialog == null) {
            controlDialog = new PeakControlDialog(IGV.getMainFrame());
        }
        controlDialog.setVisible(true);
    }


    public static float getScoreThreshold() {
        return scoreThreshold;
    }

    public static void setScoreThreshold(float t) {
        scoreThreshold = t;
        clearFilteredLists();
    }

    public static ColorOption getColorOption() {
        return colorOption;
    }

    public static void setShadeOption(ColorOption colorOption) {
        PeakTrack.colorOption = colorOption;
    }

    public static float getFoldChangeThreshold() {
        return foldChangeThreshold;
    }

    public static void setFoldChangeThreshold(float foldChangeThreshold) {
        PeakTrack.foldChangeThreshold = foldChangeThreshold;
        clearFilteredLists();
    }


    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, ReferenceFrame frame) {
        double regionScore = 0;
        int interval = end - start;
        if (interval == 0) {
            return Float.MIN_VALUE;
        }

        List<Peak> scores = getFilteredPeaks(chr);
        int startIdx = FeatureUtils.getIndexBefore(start, scores);


        for (int i = startIdx; i < scores.size(); i++) {
            Peak score = scores.get(i);
            if (score.getEnd() < start) continue;
            if (score.getStart() > end) break;
            int w = Math.min(end, score.getEnd()) - Math.max(start, score.getStart());
            float value = score.getScore();
            regionScore += value * w;
        }
        return (float) (regionScore / interval);
    }

    public Peak getPeakInstersecting(String chr, double position) {
        List<Peak> scores = getFilteredPeaks(chr);
        int startIdx = FeatureUtils.getIndexBefore(position, scores);


        for (int i = startIdx; i < scores.size(); i++) {
            Peak score = scores.get(i);
            if (score.getEnd() < position) continue;
            if (score.getStart() > position) break;
            return score;
        }
        return null;

    }


    public static boolean isShowPeaks() {
        return showPeaks;
    }

    public static void setShowPeaks(boolean b) {
        showPeaks = b;
    }

    public static boolean isShowSignals() {
        return showSignals;
    }

    public static void setShowSignals(boolean b) {
        showSignals = b;
    }

    enum ColorOption {
        SCORE, FOLD_CHANGE
    }


}
