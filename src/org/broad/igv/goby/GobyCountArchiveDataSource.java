/*
 * The MIT License (MIT)
 *  Copyright (c) 2007-2015 by Institute for Computational Biomedicine,
 *                                          Weill Medical College of Cornell University.
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

package org.broad.igv.goby;

import edu.cornell.med.icb.goby.counts.CachingCountsArchiveReader;
import edu.cornell.med.icb.goby.counts.CountBinningAdapterI;
import edu.cornell.med.icb.goby.counts.CountBinningAdaptor;
import edu.cornell.med.icb.goby.counts.CountsReader;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.data.BasicScore;
import org.broad.igv.data.CoverageDataSource;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.util.ResourceLocator;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Vector;

/**
 * A data source for <a href="http://goby.campagnelab.org">Goby</a> compressed base-level histograms (.counts files).
 *
 *
 * @author Fabien Campagne
 *         Date: 6/10/11
 *         Time: 3:18 PM
 */
public class GobyCountArchiveDataSource implements CoverageDataSource {
    static final Logger LOG = Logger.getLogger(TDFReader.class);

    CachingCountsArchiveReader counts;
    private double currentMax = 4;
    private String filename;
    private boolean someIdsStartWithChr;
    private ObjectSet<String> ids;
    private WindowFunction selectedWindowFunction;
    long numBasesSeen;
    long numSitesSeen;
    private boolean hasPrecomputedStats;
    private boolean doNormalize;
    private double normalizationFactor;

    public GobyCountArchiveDataSource(ResourceLocator locator) {
        init(locator.getPath());
    }

    static List<WindowFunction> availableFunctions = Arrays.asList(WindowFunction.mean, WindowFunction.max);

    private void init(String filename) {
        try {
            this.filename = FilenameUtils.removeExtension(filename);
            counts = new CachingCountsArchiveReader(this.filename);

            ids = counts.getIdentifiers();
            for (String id : ids) {
                if (id.startsWith("chr")) {
                    someIdsStartWithChr = true;
                    break;
                }
            }
            if (counts.isStatsParsed()) {
                this.numBasesSeen = counts.getTotalBasesSeen();
                this.numSitesSeen = counts.getTotalSitesSeen();
                hasPrecomputedStats = true;
            }
            boolean normalizeCounts = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.NORMALIZE_COVERAGE);
            setNormalize(normalizeCounts);
        } catch (IOException ex) {
            LOG.error("Error loading file: " + filename, ex);
            throw new DataLoadException("Error loading goby counts archive file: " + ex.toString(), filename);
        }
    }

    public GobyCountArchiveDataSource(File file) {
        init(file.getPath());
    }

    public double getDataMax() {
        //    LOG.info("Call in getDataMax currentMax= " + currentMax);
        return currentMax;

    }

    public double getDataMin() {
        //  LOG.info("Call in getDataMin ");
        boolean normalizeCounts = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.NORMALIZE_COVERAGE);
        setNormalize(normalizeCounts);
        return 0.0;
    }

    public List<LocusScore> getSummaryScoresForRange(String chr, int startLocation, int endLocation, int zoom) {
        if ("All".equals(chr)) return new Vector<LocusScore>();
        try {
            currentMax = 0;
            CountsReader reader = getCountsReader(chr);
            int initialStartPosition = startLocation;
            int binSize = getBinSize(startLocation, endLocation);

            /*       LOG.info(String.format("Call in getSummaryScoresForRange %d-%d zoom %d binSize=%d ", startLocation, endLocation,
                       zoom, binSize));
            */
            if (reader == null) {
                return null;
            }
            updateNormalizationFactor();
            CountBinningAdapterI binAdaptor = new CountBinningAdaptor(reader, binSize);
            if (counts.hasIndex()) {
                // we can only reposition if the countsreader has an index. Otherwise, we trust the
               // caching counts  archive  to have created a new reader positioned at the start of the counts sequence.
                binAdaptor.reposition(startLocation);
            }
            binAdaptor.skipTo(startLocation);

            int position = binAdaptor.getPosition();
            ObjectArrayList<LocusScore> result = new ObjectArrayList<LocusScore>();
            result.add(new BasicScore(initialStartPosition, position, 0.0f));
            while (binAdaptor.hasNextTransition() && position < endLocation) {
                binAdaptor.nextTransition();
                // position is the zero-based position before the count is changed by the transition
                position = binAdaptor.getPosition();
                // count is how many reads cover the length bases that follow position.
                final double count = selectedWindowFunction == WindowFunction.mean ? binAdaptor.getAverage() : binAdaptor.getMax();

                final double normalizedCount = count / normalizationFactor;
                final int length = binAdaptor.getLength();
                //      System.out.printf("adding results %d-%d count=%g %n", position, position + length , normalizedCount);

                BasicScore bs = new BasicScore(position, position + length, (float) normalizedCount);
                result.add(bs);

                this.currentMax = Math.max(currentMax, normalizedCount);
                if (!hasPrecomputedStats) {
                    // estimate on the fly from the subset of data seen so far:
                    numBasesSeen += length * count;
                    if (count != 0) {
                        numSitesSeen += length;
                    }
                }
            }
            result.add(new BasicScore(position, endLocation, 0.0f));

            return result;
        } catch (IOException e) {
            LOG.error(e);
            throw new DataLoadException(
                    String.format("Error getting summary scores for range %s:%d-%d in goby counts archive file %s %n",
                            chr, startLocation, endLocation, filename), filename);

        }
    }

    private int getBinSize(int startLocation, int endLocation) {
        final int regionLength = endLocation - startLocation;
        return Math.max(1, ((regionLength / 2000)));
    }

    private void updateNormalizationFactor() {
        if (numSitesSeen == 0) {
            normalizationFactor = 1.0;
        } else {
            // normalization factor is estimated such that the average value will be ~2 for autosomes, corresponding
            // to two copies.
            normalizationFactor = doNormalize ?  ((double)numBasesSeen / (double) numSitesSeen) /2.0:
                    1.0;
        }
        //  System.out.printf("normalization factor=%g%n", normalizationFactor);

    }

    private CountsReader getCountsReader(String chr) throws IOException {
        if (ids.contains(chr)) {
            return counts.getCountReader(chr);
        }
        if (!someIdsStartWithChr) {

            CountsReader reader = counts.getCountReader(chr.replaceFirst("chr", ""));
            if (reader != null) return reader;
        }
        return counts.getCountReader(chr);
    }

    public TrackType getTrackType() {
        return TrackType.COVERAGE;
    }

    public void setWindowFunction(WindowFunction statType) {
        selectedWindowFunction = statType;
        //  LOG.info("Call in setWindowFunction ");
    }

    public boolean isLogNormalized() {
        //  LOG.info("Call in isLogNormalized ");
        return false;
    }

    public void refreshData(long timestamp) {
        //     LOG.info("Call in refreshData ");
    }

    public WindowFunction getWindowFunction() {
        //     LOG.info("Call in getWindowFunction ");
        return selectedWindowFunction;
    }

    public Collection<WindowFunction> getAvailableWindowFunctions() {
        //  LOG.info("Call in getAvailableWindowFunctions ");
        return availableFunctions;
    }

    @Override
    public void dispose() {

    }

    public String getPath() {
        return filename;
    }

    public void setNormalize(boolean normalize) {
        this.doNormalize = normalize;
    }

    public boolean getNormalize(){
        return this.doNormalize;
    }
}
