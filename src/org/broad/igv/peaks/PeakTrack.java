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

package org.broad.igv.peaks;

import org.broad.igv.data.DataSource;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.tdf.TDFDataSource;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.NamedRunnable;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.io.BufferedReader;
import java.io.IOException;
import java.lang.ref.SoftReference;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 * @date Apr 22, 2011
 */
public class PeakTrack extends AbstractTrack {

    private static org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(PeakTrack.class);

    static List<SoftReference<PeakTrack>> instances = new ArrayList();

    private static PeakControlDialog controlDialog;
    private static float scoreThreshold = 30;
    private static float foldChangeThreshold = 0;

    private static ColorOption colorOption = ColorOption.SCORE;
    private static boolean showPeaks = true;
    private static boolean showSignals = false;
    static boolean commandBarAdded = false;
    static int timeStep = 0;
    static boolean animate = false;

    int nTimePoints;


    Map<String, List<Peak>> peakMap = new HashMap();
    Map<String, List<Peak>> filteredPeakMap = new HashMap();
    Renderer renderer = new PeakRenderer();


    // Path to the compressed signal (TDF) file and data source
    String signalPath;
    private WrappedDataSource signalSource;

    // Paths to the time series signal files and data sources
    String[] timeSignalPaths;
    WrappedDataSource[] timeSignalSources;

    // Data range
    DataRange scoreDataRange = new DataRange(0, 0, 100);
    DataRange signalDataRange = new DataRange(0, 0, 1000f);

    int bandHeight;
    int signalHeight;
    int peakHeight;
    int gapHeight;

    Genome genome;
    private String peaksPath;

    PeakParser parser;

    static synchronized boolean isCommandBarAdded() {
        boolean retValue = commandBarAdded;
        commandBarAdded = true;
        return retValue;
    }


    /**
     * @param locator -- path to a peaks.cfg file
     * @param genome
     * @throws IOException
     */
    public PeakTrack(final ResourceLocator locator, Genome genome) throws IOException {
        super(locator);
        this.genome = genome;
        setHeight(20);


        try {
            long t0 = System.currentTimeMillis();
            parser = new PeakParser(locator.getPath());
            this.getAllPeaks("chr2");
            long dt = System.currentTimeMillis() - t0;
            //log.info("Loaded bin: " + locator.getPath() + ": " + dt);

            TrackProperties props = new TrackProperties();
            ParsingUtils.parseTrackLine(parser.trackLine, props);
            setProperties(props);

            nTimePoints = parser.nTimePoints;
            signalPath = parser.signalsPath;
            timeSignalPaths = parser.timeSignalsPath;


        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }


        instances.add(new SoftReference(this));

        if (!isCommandBarAdded()) {
            IGV.getInstance().getContentPane().addCommandBar(new PeakCommandBar());

        }
    }


    /**
     * timePoints=0,30,60,120
     * peaks=http://www.broadinstitute.org/igvdata/ichip/peaks/AHR.peak
     * signals=http://www.broadinstitute.org/igvdata/ichip/tdf/compressed/AHR.merged.bam.tdf
     * timeSignals=http://www.broadinstitute.org/igvdata/ichip/tdf/timecourses/AHR_0/AHR_0.merged.bam.tdf,http...
     *
     * @param path
     * @throws IOException
     */

    private void parse(String path) throws IOException {

        BufferedReader br = null;


        try {
            br = ParsingUtils.openBufferedReader(path);

            String nextLine = br.readLine();
            if (nextLine.startsWith("track")) {
                TrackProperties props = new TrackProperties();
                ParsingUtils.parseTrackLine(nextLine, props);
                setProperties(props);
            }

            nextLine = br.readLine();
            String[] tokens = nextLine.split("=");
            if (tokens.length < 2 || !tokens[0].equals("timePoints")) {
                throw new RuntimeException("Unexpected timePoints line: " + nextLine);
            }
            tokens = tokens[1].split(",");
            nTimePoints = tokens.length;

            nextLine = br.readLine();
            tokens = nextLine.split("=");
            if (tokens.length < 2 || !tokens[0].equals("peaks")) {
                throw new RuntimeException("Unexpected timePoints line: " + nextLine);
            }
            peaksPath = tokens[1];


            nextLine = br.readLine();
            tokens = nextLine.split("=");
            if (tokens.length < 2 || !tokens[0].equals("signals")) {
                throw new RuntimeException("Unexpected timePoints line: " + nextLine);
            }
            signalPath = tokens[1];

            nextLine = br.readLine();
            tokens = nextLine.split("=");
            if (tokens.length < 2 || !tokens[0].equals("timeSignals")) {
                throw new RuntimeException("Unexpected timePoints line: " + nextLine);
            }
            timeSignalPaths = tokens[1].split(",");


        } finally {
            if (br != null) br.close();
        }


    }


    @Override
    public IGVPopupMenu getPopupMenu(TrackClickEvent te) {
        return new PeakTrackMenu(this, te);
    }

    @Override
    public DataRange getDataRange() {
        return showSignals ? signalDataRange : scoreDataRange;
    }

    @Override
    public void setDataRange(DataRange axisDefinition) {
        if (showSignals) {
            signalDataRange = axisDefinition;
        } else {
            scoreDataRange = axisDefinition;
        }
    }

    @Override
    public void load(RenderContext context) {
        try {
            getFilteredPeaks(context.getChr());
        } catch (IOException e) {
            log.error("Error loading peaks", e);
        }
    }


    public void render(RenderContext context, Rectangle rect) {

        try {
            List<Peak> peakList = getFilteredPeaks(context.getChr());
            if (peakList == null) {
                return;
            }

            renderer.render(peakList, context, rect, this);

        } catch (IOException e) {
            log.error("Error loading peaks", e);
        }
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
            return nTimePoints * h + gapHeight;
        }
    }

    @Override
    public void setHeight(int h, boolean force) {
        super.setHeight(h, force);

        int nBands = getDisplayMode() == DisplayMode.COLLAPSED ? 1 : nTimePoints;

        bandHeight = h / nBands;
        peakHeight = Math.max(5, Math.min(bandHeight / 3, 10));
        signalHeight = bandHeight - peakHeight - gapHeight;

    }

    @Override
    public void setDisplayMode(DisplayMode mode) {
        super.setDisplayMode(mode);
        if (mode == Track.DisplayMode.EXPANDED) {
            setHeight(nTimePoints * bandHeight + gapHeight);
        } else {
            setHeight(bandHeight);
        }
    }


    public static void setAnimate(boolean animate) {
        PeakTrack.animate = animate;
        if(animate) {
            startAnimationThread();
        }
    }

    private static synchronized void startAnimationThread() {

        timeStep = 0;

        Runnable runnable = new Runnable() {
            @Override
            public void run() {

                while (animate) {
                    try {
                        Thread.sleep(1000) ;
                        timeStep ++;
                        if(timeStep == 4) {
                            Thread.sleep(500);
                            timeStep = 0;
                        }
                        IGV.getInstance().doRefresh();

                    } catch (InterruptedException e) {
                        setAnimate(false);
                        return;
                    }
                }
            }
        };

        (new Thread(runnable)).start();

    }


    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {
        try {
            boolean foundValue = false;
            StringBuffer buf = new StringBuffer();
            buf.append(getName());
            buf.append("<br>");
            if (showPeaks) {
                List<Peak> scores = getFilteredPeaks(chr);
                LocusScore score = getLocusScoreAt(scores, position, frame);
                if (score != null) {
                    foundValue = true;
                    buf.append(score.getValueString(position, getWindowFunction()));
                    if (showSignals) {
                        buf.append("<br>");
                    }
                }
            }

            final WrappedDataSource signalSource = getSignalSource();
            if (showSignals && signalSource != null) {
                List<LocusScore> scores = signalSource.getSummaryScoresForRange(chr, (int) frame.getOrigin(), (int) frame.getEnd(), frame.getZoom());
                LocusScore score = getLocusScoreAt(scores, position, frame);
                if (score != null) {
                    foundValue = true;
                    buf.append("Score = " + score.getScore());
                }
            }
            return foundValue ? buf.toString() : null;
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            return "Error loading peaks: " + e.toString();
        }
    }


    // TODO -- the code below is an exact copy of code in DataTrack.   Refactor to share this.

    private LocusScore getLocusScoreAt(List<? extends LocusScore> scores, double position, ReferenceFrame frame) {

        if (scores == null) {
            return null;
        } else {
            // give a 2 pixel window, otherwise very narrow features will be missed.
            double bpPerPixel = frame.getScale();
            int buffer = (int) (2 * bpPerPixel);    /* * */
            return FeatureUtils.getFeatureAt(position, buffer, scores);
        }
    }

    public synchronized List<Peak> getFilteredPeaks(String chr) throws IOException {
        List<Peak> filteredPeaks = filteredPeakMap.get(chr);
        if (filteredPeaks == null) {
            filteredPeaks = new ArrayList();
            List<Peak> allPeaks = getAllPeaks(chr);
            if (allPeaks != null) {
                for (Peak peak : allPeaks) {
                    if (peak.getCombinedScore() >= scoreThreshold &&
                            peak.getFoldChange() >= foldChangeThreshold) {
                        filteredPeaks.add(peak);
                    }
                }
            }
        }
        filteredPeakMap.put(chr, filteredPeaks);


        return filteredPeaks;
    }

    private List<Peak> getAllPeaks(String chr) throws IOException {

        List<Peak> peaks = peakMap.get(chr);
        if (peaks == null) {
            peaks = parser.loadPeaks(chr);
            peakMap.put(chr, peaks);
        }
        return peaks;
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




    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, String frameName) {

        int interval = end - start;
        if (interval <= 0) {
            return Float.MIN_VALUE;
        }

        try {
            List<Peak> scores = getFilteredPeaks(chr);
            int startIdx = Math.max(0, FeatureUtils.getIndexBefore(start, scores));

            float regionScore = Float.MIN_VALUE;
            for (int i = startIdx; i < scores.size(); i++) {
                Peak score = scores.get(i);
                if (score.getEnd() < start) continue;
                if (score.getStart() > end) break;
                final float v = score.getScore();
                if (v > regionScore) regionScore = v;
            }
            return regionScore;
        } catch (IOException e) {
            return Float.MIN_VALUE;
        }
    }


    /**
     * Get the closet filter peak, within 2kb, of the given position.
     *
     * @param chr
     * @param position
     * @return
     */
    public Peak getFilteredPeakNearest(String chr, double position) {
        try {
            List<Peak> scores = getFilteredPeaks(chr);
            int startIdx = FeatureUtils.getIndexBefore(position, scores);

            Peak closestPeak = null;
            double closestDistance = Integer.MAX_VALUE;
            if (startIdx >= 0) {
                if (startIdx > 0) startIdx--;
                for (int i = startIdx; i < scores.size(); i++) {
                    Peak peak = scores.get(i);
                    if (position > peak.getStart() && position < peak.getEnd()) {
                        return peak;
                    }
                    double distance = Math.min(Math.abs(position - peak.getStart()), Math.abs(position - peak.getEnd()));
                    if (distance > closestDistance) {
                        return closestDistance < 2000 ? closestPeak : null;

                    } else {
                        closestDistance = distance;
                        closestPeak = peak;
                    }

                }
            }
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
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

    public int getTimeStep() {
        return timeStep;
    }

    public DataSource[] getTimeSignalSources() {

        if (timeSignalSources == null) {
            if (timeSignalPaths != null && timeSignalPaths.length > 0) {
                timeSignalSources = new WrappedDataSource[timeSignalPaths.length];
                for (int i = 0; i < timeSignalPaths.length; i++) {
                    try {
                        timeSignalSources[i] = new WrappedDataSource(new TDFDataSource(TDFReader.getReader(timeSignalPaths[i]), 0, "", genome));
                        timeSignalSources[i].setNormalizeCounts(true, 1.0e9f);
                    } catch (Exception e) {
                        timeSignalSources[i] = null;
                        e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                    }
                }
            }
        }

        return timeSignalSources;
    }

    boolean signalSourceLoading = false;


    public WrappedDataSource getSignalSource() {
        return signalSource;
    }

    /**
     * Called by Renderer.  Loads signal source, if needed, and forces load of specified data interval.  This is
     * a bit of a hack.
     *
     * @param chr
     * @param contextStart
     * @param contextEnd
     * @param zoom
     * @return
     */
    WrappedDataSource getSignalSource(final String chr, final int contextStart, final int contextEnd, final int zoom) {
        if (signalSource == null && signalPath != null && !signalSourceLoading) {
            signalSourceLoading = true;
            NamedRunnable runnable = new NamedRunnable() {
                public void run() {
                    signalSource = new WrappedDataSource(new TDFDataSource(TDFReader.getReader(signalPath), 0, "", genome));
                    signalSource.setNormalizeCounts(true, 1.0e9f);
                    signalSource.getSummaryScoresForRange(chr, contextStart, contextEnd, zoom);
                }

                public String getName() {
                    return "Load " + signalPath;
                }
            };
            LongRunningTask.submit(runnable);
        }
        return signalSource;
    }


    public String getSignalPath() {
        return signalPath;
    }

    enum ColorOption {
        SCORE, FOLD_CHANGE
    }


    class WrappedDataSource implements DataSource {

        TDFDataSource source;

        WrappedDataSource(TDFDataSource source) {
            this.source = source;
        }

        public List<LocusScore> getSummaryScoresForRange(String chr, int startLocation, int endLocation, int zoom) {

            List<LocusScore> scores = new ArrayList(1000);


            if (scoreThreshold <= 0 && foldChangeThreshold <= 0) {
                return source.getSummaryScoresForRange(chr, startLocation, endLocation, zoom);
            } else {
                try {
                    List<Peak> peaks = getFilteredPeaks(chr);
                    if (peaks == null) {
                        return scores;
                    }
                    int startIdx = FeatureUtils.getIndexBefore(startLocation, peaks);
                    if (startIdx >= 0) {
                        for (int i = startIdx; i < peaks.size(); i++) {
                            Peak peak = peaks.get(i);

                            final int peakEnd = peak.getEnd();
                            if (peakEnd < startLocation) continue;

                            final int peakStart = peak.getStart();
                            if (peakStart > endLocation) break;

                            List<LocusScore> peakScores = source.getSummaryScoresForRange(chr, peakStart, peakEnd, zoom);
                            for (LocusScore ps : peakScores) {
                                if (ps.getEnd() < peakStart) continue;
                                if (ps.getStart() > peakEnd) break;
                                scores.add(ps);
                            }

                        }
                    }
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
                return scores;
            }
        }


        public double getDataMax() {
            return source.getDataMax();
        }

        public double getDataMin() {
            return source.getDataMin();
        }

        public TrackType getTrackType() {
            return source.getTrackType();
        }

        public void setWindowFunction(WindowFunction statType) {
            source.setWindowFunction(statType);
        }

        public boolean isLogNormalized() {
            return source.isLogNormalized();
        }

        public WindowFunction getWindowFunction() {
            return source.getWindowFunction();
        }

        public Collection<WindowFunction> getAvailableWindowFunctions() {
            return source.getAvailableWindowFunctions();
        }

        @Override
        public void dispose() {
            source.dispose();
        }

        public void setNormalizeCounts(boolean b, float v) {
            source.setNormalizeCounts(b, v);
        }

        public void updateGenome(Genome genome) {
            source.updateGenome(genome);
        }
    }


    private InViewInterval computeScale(double origin, double end, List<LocusScore> scores) {

        InViewInterval interval = new InViewInterval();

        if (scores.size() == 1) {
            interval.dataMax = Math.max(0, scores.get(0).getScore());
            interval.dataMin = Math.min(0, scores.get(0).getScore());
        } else {
            interval.startIdx = 0;
            interval.endIdx = scores.size();
            for (int i = 1; i < scores.size(); i++) {
                if (scores.get(i).getEnd() >= origin) {
                    interval.startIdx = i - 1;
                    break;
                }
            }

            for (int i = interval.startIdx + 1; i < scores.size(); i++) {
                LocusScore locusScore = scores.get(i);
                float value = locusScore.getScore();
                if (Float.isNaN(value)) value = 0;
                interval.dataMax = Math.max(interval.dataMax, value);
                interval.dataMin = Math.min(interval.dataMin, value);
                if (locusScore.getStart() > end) {
                    interval.endIdx = i;
                    break;
                }
            }
        }

        return interval;
    }

    class InViewInterval {
        int startIdx;
        int endIdx;
        float dataMax = 0;
        float dataMin = 0;
    }


}
