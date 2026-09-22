package org.igv.alignment;

import org.apache.commons.math3.stat.StatUtils;
import org.igv.logging.*;
import org.igv.util.collections.DownsampledDoubleArrayList;

import java.util.Arrays;

/**
 * Insert size ("TLEN") and pair orientation statistics for a single library.
 * <p>
 * Samples are pooled over every load, and no single load can contribute more than a small fraction of the
 * pool, so the thresholds do not depend on which region happened to be loaded first.  Loads run in
 * parallel, so each accumulates into its own instance (see {@link #forLoad}) and merges that into the
 * shared instance for the library when it completes, rather than updating shared state read by read.
 *
 * @author jrobinso
 * @date Mar 11, 2011
 */
public class PEStats {

    private static Logger log = LogManager.getLogger(PEStats.class);

    public enum Orientation {FR, RF, F1F2, F2F1}

    String library;


    //Maximum number of insertSizes to store
    private static final int MAX = 20000;

    // Maximum number of insertSizes contributed by a single load.  A deep load can supply more proper pairs
    // than the whole pool holds, so without this cap a single region could set the thresholds for the
    // entire session, which is the behavior this class exists to avoid.
    private static final int MAX_PER_LOAD = 1000;

    // Minimum number of insertSizes required before thresholds are computed
    private static final int MIN_SAMPLES = 100;

    // Trim at this many robust standard deviations before taking a percentile.  Wide enough that a well
    // behaved library loses nothing, narrow enough to separate out a discordant population.
    private static final int TRIM_DEVIATIONS = 10;

    private DownsampledDoubleArrayList insertSizes;
    private volatile int minThreshold = 10;
    private volatile int maxThreshold = 5000;

    // Orientation counts
    int frCount = 0;
    int rfCount = 0;

    int f1f2Count = 0;
    int f2f1Count = 0;

    int totalCount = 0;

    volatile Orientation orientation = Orientation.FR;


    /**
     * For paired arc view which are outside of the midrange
     * TODO Allow user to set?
     */
    private static final int minOutlierInsertSizePercentile = 95;
    private static final int maxOutlierInsertSizePercentile = 5;

    private volatile int minOutlierInsertSize = minThreshold;
    private volatile int maxOutlierInsertSize = maxThreshold;

    public PEStats(String library) {
        this(library, MAX);
    }

    /**
     * Return an instance for accumulating the sample of a single load.  Its sample is capped at
     * MAX_PER_LOAD and merged into the shared instance for the library when the load completes.
     */
    static PEStats forLoad(String library) {
        return new PEStats(library, MAX_PER_LOAD);
    }

    private PEStats(String library, int maxSamples) {
        this.library = library;
        this.insertSizes = new DownsampledDoubleArrayList(100, maxSamples);
    }


    public synchronized void update(Alignment alignment) {

        if (alignment.isProperPair()) {
            insertSizes.add(Math.abs(alignment.getInferredInsertSize()));
            String po = alignment.getPairOrientation();
            if (po != null && po.length() == 4) {
                if (po.charAt(0) == 'F') {
                    if (po.charAt(2) == 'F') {
                        if (po.charAt(1) == '1') {
                            f1f2Count++;
                        } else {
                            f2f1Count++;
                        }
                    } else if (po.charAt(2) == 'R') {
                        frCount++;

                    }
                } else if (po.charAt(0) == 'R') {
                    if (po.charAt(2) == 'F') {
                        rfCount++;
                    } else if (po.charAt(2) == 'R') {
                        if (po.charAt(1) == '1') {
                            f2f1Count++;
                        } else {
                            f1f2Count++;
                        }
                    }
                }
            }
            totalCount++;
        }
    }

    /**
     * Merge the statistics accumulated by a single load.  The merged instance is not shared with any
     * other thread, so only this instance needs guarding.
     *
     * @param other
     */
    public synchronized void merge(PEStats other) {

        for (double isize : other.insertSizes.toArray()) {
            insertSizes.add(isize);
        }
        frCount += other.frCount;
        rfCount += other.rfCount;
        f1f2Count += other.f1f2Count;
        f2f1Count += other.f2f1Count;
        totalCount += other.totalCount;
    }

    public synchronized void computeInsertSize(double minPercentile, double maxPercentile) {

        // A minimum percentile of zero means no read is to be colored as a small insert.  It is not a
        // percentile an estimate can be taken at, and it applies before there is a sample to estimate from.
        final boolean noMinimum = minPercentile <= 0;
        if (noMinimum) {
            minThreshold = 0;
        }

        if (insertSizes.size() > MIN_SAMPLES) {

            final double[] sample = trimOutliers(insertSizes.toArray());

            if (!noMinimum) {
                minThreshold = computePercentile(sample, minPercentile);
            }
            maxThreshold = computePercentile(sample, maxPercentile);

            minOutlierInsertSize = computePercentile(sample, minOutlierInsertSizePercentile);
            maxOutlierInsertSize = computePercentile(sample, maxOutlierInsertSizePercentile);
        }
    }

    public int getMinThreshold() {
        return minThreshold;
    }

    public int getMaxThreshold() {
        return maxThreshold;
    }

    public Orientation getOrientation() {
        if (orientation == null) {
            computeExpectedOrientation();
        }
        return orientation;
    }

    public synchronized void computeExpectedOrientation() {

        if(totalCount > 100) {
            int ffCount = f1f2Count + f2f1Count;
            if (ffCount > frCount && ffCount > rfCount) {
                if (f1f2Count > f2f1Count) {
                    orientation = Orientation.F1F2;
                } else {
                    orientation = Orientation.F2F1;
                }
            } else if (rfCount > frCount && rfCount > ffCount) {
                orientation = Orientation.RF;
            } else {
                orientation = Orientation.FR;
            }
        }
        else {
            orientation = Orientation.FR;
        }
    }

    private static int computePercentile(double[] sample, double percentile) {
        return (int) StatUtils.percentile(sample, 0, sample.length, percentile);
    }

    /**
     * Return the sample with its outlier population, if any, removed.  The percentiles used for the
     * thresholds sit in the tails of the distribution, so a small fraction of discordant pairs -- which
     * form a separate population far from the mode -- is enough to drag a threshold past every read worth
     * flagging.  That population is located with the median and the median absolute deviation, neither of
     * which a minority of outliers can shift.  A library with no such population loses nothing, so
     * ordinary data is unaffected.
     *
     * @param insertSizes
     */
    private static double[] trimOutliers(double[] insertSizes) {

        double[] sorted = insertSizes.clone();
        Arrays.sort(sorted);

        double median = medianOfSorted(sorted);
        double[] deviations = new double[sorted.length];
        for (int i = 0; i < sorted.length; i++) {
            deviations[i] = Math.abs(sorted[i] - median);
        }
        Arrays.sort(deviations);
        double scale = 1.4826 * medianOfSorted(deviations);
        if (scale <= 0) {
            return sorted;
        }

        double low = median - TRIM_DEVIATIONS * scale;
        double high = median + TRIM_DEVIATIONS * scale;
        int from = 0;
        while (from < sorted.length && sorted[from] < low) from++;
        int to = sorted.length;
        while (to > from && sorted[to - 1] > high) to--;

        // Refuse to trim away a majority.  If that happens the sample is not a contaminated unimodal
        // distribution, and there is no basis for calling any part of it an outlier.
        return 2 * (to - from) >= sorted.length ? Arrays.copyOfRange(sorted, from, to) : sorted;
    }

    private static double medianOfSorted(double[] sorted) {
        int n = sorted.length;
        return n % 2 == 1 ? sorted[n / 2] : (sorted[n / 2 - 1] + sorted[n / 2]) / 2;
    }

    int getMinOutlierInsertSize() {
        return minOutlierInsertSize;
    }

    int getMaxOutlierInsertSize() {
        return maxOutlierInsertSize;
    }
}
