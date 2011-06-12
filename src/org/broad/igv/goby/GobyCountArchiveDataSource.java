package org.broad.igv.goby;

import edu.cornell.med.icb.goby.counts.*;
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
 * @author Fabien Campagne
 *         Date: 6/10/11
 *         Time: 3:18 PM
 */
public class GobyCountArchiveDataSource implements CoverageDataSource {
    static final Logger LOG = Logger.getLogger(TDFReader.class);

    CountsArchiveReader counts;
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
            counts = new CountsArchiveReader(this.filename);
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
  //      LOG.info("Call in getDataMin ");
        boolean normalizeCounts = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.NORMALIZE_COVERAGE);
        setNormalize(normalizeCounts);
        return 0.0;
    }

    public List<LocusScore> getSummaryScoresForRange(String chr, int startLocation, int endLocation, int zoom) {
        if ("All".equals(chr)) return new Vector<LocusScore>();
        try {

            CountsReader reader = getCountsReader(chr);

            final int regionLength = endLocation - startLocation;
            int binSize = Math.max(1,((regionLength / 2000)));

        /*    LOG.info(String.format("Call in getSummaryScoresForRange %d-%d zoom %d binSize=%d ", startLocation, endLocation,
                    zoom, binSize));
          */
            if (reader == null) {
                return new Vector<LocusScore>();
            }
            updateNormalizationFactor();

          //  final boolean adapterChoice = regionLength > 50000 && regionLength < 1000000;
            final boolean adapterChoice = regionLength > 3000000 && regionLength < 10000000;
      //      System.out.println(adapterChoice?"CountAdaptiveBinningAdaptor":"CountBinningAdaptor");
            CountBinningAdapterI binAdaptor = adapterChoice ?new CountAdaptiveBinningAdaptor(reader):
                    new CountBinningAdaptor(reader,binSize);

            binAdaptor.skipTo(startLocation);
            ObjectArrayList<LocusScore> result = new ObjectArrayList<LocusScore>();
            int position = binAdaptor.getPosition();
            while (binAdaptor.hasNextTransition() && position < endLocation) {
                binAdaptor.nextTransition();
                // position is the zero-based position before the count is changed by the transition
                position = binAdaptor.getPosition();
                // count is how many reads cover the length bases that follow position.
                final double count = selectedWindowFunction == WindowFunction.mean ? binAdaptor.getAverage() : binAdaptor.getMax();

                final double normalizedCount = count / normalizationFactor;
                final int length = binAdaptor.getLength();
          //      System.out.printf("adding results %d-%d count=%g %n", position, position + length , normalizedCount);
                BasicScore bs = new BasicScore(position, position + length , (float) normalizedCount);

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

         //   System.out.printf("returning %d results %n", result.size());
            return result;
        } catch (IOException e) {
            LOG.error(e);
            throw new DataLoadException(
                    String.format("Error getting summary scores for range %s:%d-%d in goby counts archive file %s %n",
                            chr, startLocation, endLocation, filename), filename);

        }
    }

    private void updateNormalizationFactor() {
        if (numSitesSeen == 0) {
            normalizationFactor = 1.0;
        } else {
            normalizationFactor = doNormalize ? (double) numBasesSeen / (double) numSitesSeen :
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

    public String getPath() {
        return filename;
    }

    public void setNormalize(boolean normalize) {
        this.doNormalize = normalize;
    }
}
