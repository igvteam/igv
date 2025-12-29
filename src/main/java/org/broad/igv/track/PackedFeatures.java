package org.broad.igv.track;

import org.broad.igv.logging.*;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.Strand;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import htsjdk.tribble.Feature;

import java.util.*;

import static org.broad.igv.feature.FeatureUtils.FEATURE_CENTER_COMPARATOR;

/**
 * Represents a table of features, packed so there is no overlap.
 * Features are packed into rows, accessible via {@link #getRows}
 *
 * @author jrobinso
 * @date Oct 7, 2010
 */
public class PackedFeatures<T extends Feature> {
    protected String trackName;
    protected String chr;
    protected int start;
    protected int end;
    protected List<T> features;
    protected List<FeatureRow> rows;
    private static Logger log = LogManager.getLogger(PackedFeatures.class);
    protected int maxFeatureLength = 0;
    protected static int maxLevels = 1000000;
    private Track.DisplayMode displayMode;
    private boolean groupByStrand;

    /**
     * Cache for storing features sorted by center position.   This is used for hotkey feature jumping.
     */
    private List<T> centerSortedFeatures;


    /**
     * Create an empty PackedFeatures object.
     *
     * @param chr
     * @param start
     * @param end
     */
    PackedFeatures(String chr, int start, int end) {
        this.chr = chr;
        this.start = start;
        this.end = end;
        features = Collections.emptyList();
        rows = Collections.emptyList();
    }

    PackedFeatures(String chr, int start, int end, Iterator<T> iter, Track.DisplayMode displayMode, boolean groupByStrand) {
        this.chr = chr;
        this.start = start;
        this.end = end;
        features = new ArrayList<>(100);
        while (iter.hasNext()) {
            T feature = (T) iter.next();
            maxFeatureLength = Math.max(maxFeatureLength,
                    getFeatureEndForPacking(feature) - getFeatureStartForPacking(feature));
            features.add(feature);
        }

        pack(displayMode, groupByStrand);
    }

    public void pack(Track.DisplayMode displayMode, boolean groupByStrand) {

        // Repack if groupByStrand has changed, or if display mode has switched to/from COLLAPSED
        if (this.displayMode == null ||
                this.groupByStrand != groupByStrand ||
                (this.displayMode == Track.DisplayMode.COLLAPSED && displayMode != Track.DisplayMode.COLLAPSED) ||
                (this.displayMode != Track.DisplayMode.COLLAPSED && displayMode == Track.DisplayMode.COLLAPSED)) {

            this.displayMode = displayMode;
            this.groupByStrand = groupByStrand;


            if (groupByStrand) {
                List<T> posFeatures = new ArrayList<>();
                List<T> negFeatures = new ArrayList<>();
                for (T f : features) {
                    if (f instanceof IGVFeature && ((IGVFeature) f).getStrand() == Strand.NEGATIVE) {
                        negFeatures.add(f);
                    } else {
                        posFeatures.add(f);
                    }
                }
                if (displayMode == Track.DisplayMode.COLLAPSED) {
                    rows = Arrays.asList(new FeatureRow(posFeatures), new FeatureRow(negFeatures));
                } else {
                    rows = packFeatures(posFeatures.iterator());
                    rows.addAll(packFeatures(negFeatures.iterator()));
                }
            } else {

                if (displayMode == Track.DisplayMode.COLLAPSED) {
                    rows = Arrays.asList(new FeatureRow(features));
                } else {
                    rows = packFeatures(features.iterator());
                }
            }
        }
    }


    /**
     * Some types of Features (splice junctions) should be packed on the same row even if start and end overlap.
     * This can be overridden in a subclass
     *
     * @param feature
     * @return
     */
    protected int getFeatureStartForPacking(Feature feature) {
        return feature.getStart();
    }


    /**
     * Some types of Features (splice junctions) should be packed on the same row even if start and end overlap.
     * This can be overridden in a subclass
     *
     * @param feature
     * @return
     */
    protected int getFeatureEndForPacking(Feature feature) {
        return feature.getEnd();
    }

    int getRowCount() {
        return rows.size();
    }

    public boolean containsInterval(String chr, int start, int end) {
        return this.getChr().equals(chr) && start >= this.getStart() && end <= this.getEnd();
    }

    public boolean overlapsInterval(String chr, int start, int end) {
        return this.getChr().equals(chr) && start <= this.end && end >= this.start;
    }

    /**
     * Allocates each feature to the rows such that there is no overlap.
     *
     * @param iter TabixLineReader wrapping the collection of alignments. Note that this should
     *             really be an Iterator<T>, but it can't be subclassed if that's the case.
     */
    List<FeatureRow> packFeatures(Iterator iter) {

        List<FeatureRow> rows = new ArrayList(10);
        if (iter == null || !iter.hasNext()) {
            return rows;
        }

        int totalCount = 0;

        LinkedHashMap<Integer, PriorityQueue<T>> bucketArray = new LinkedHashMap();
        Comparator pqComparator = (Comparator<T>) (row1, row2) -> (row2.getEnd() - row2.getStart()) - (row1.getEnd() - row2.getStart());

        // Allocate features to buckets,  1 bucket per base position
        while (iter.hasNext()) {
            T feature = (T) iter.next();

            int bucketNumber = getFeatureStartForPacking(feature);

            PriorityQueue<T> bucket = bucketArray.get(bucketNumber);
            if (bucket == null) {
                bucket = new PriorityQueue<T>(5, pqComparator);
                bucketArray.put(bucketNumber, bucket);
            }
            bucket.add(feature);
            totalCount++;

        }

        // Allocate features to rows, pulling at most 1 per bucket for each row
        FeatureRow currentRow = new FeatureRow();
        int allocatedCount = 0;
        int nextStart = Integer.MIN_VALUE;

        int lastAllocatedCount = -1;
        while (allocatedCount < totalCount && rows.size() < maxLevels) {

            // Check to prevent infinite loops
            if (lastAllocatedCount == allocatedCount) {

                if (IGV.hasInstance()) {
                    String msg = "Infinite loop detected while packing features for track: " +
                            ".<br>Not all features will be shown." +
                            "<br>Please contact igv-team@broadinstitute.org";

                    log.error(msg);
                    MessageUtils.showMessage(msg);
                }
                break;
            }
            lastAllocatedCount = allocatedCount;

            // Next row Loop through alignments until we reach the end of the interval

            PriorityQueue<T> bucket = null;
            // Advance to nextLine occupied bucket

            ArrayList<Integer> emptyBucketKeys = new ArrayList();
            for (Integer key : bucketArray.keySet()) {
                if (key >= nextStart) {
                    bucket = bucketArray.get(key);

                    T feature = bucket.poll();

                    if (bucket.isEmpty()) {
                        emptyBucketKeys.add(key);
                    }
                    currentRow.addFeature(feature);
                    nextStart = currentRow.end + FeatureTrack.MINIMUM_FEATURE_SPACING;
                    allocatedCount++;
                }
            }
            for (Integer key : emptyBucketKeys) {
                bucketArray.remove(key);
            }


            // We've reached the end of the interval,  start a new row
            if (currentRow.features.size() > 0) {
                rows.add(currentRow);
                lastAllocatedCount = 0;
            }
            currentRow = new FeatureRow();
            nextStart = 0;
        }
        // Add the last row
        if (currentRow.features.size() > 0) {
            rows.add(currentRow);
        }

        return rows;
    }

    public String getChr() {
        return chr;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public List<T> getFeatures() {
        return features;
    }

    public List<FeatureRow> getRows() {
        return rows;
    }

    public int getMaxFeatureLength() {
        return maxFeatureLength;
    }

    public List<T> getCenterSortedFeatures() {
        if(centerSortedFeatures == null && features != null) {
            centerSortedFeatures = new ArrayList<>(features);
            Collections.sort(centerSortedFeatures, FEATURE_CENTER_COMPARATOR);
        }
        return centerSortedFeatures;
    }

    public void setCenterSortedFeatures(List<T> centerSortedFeatures) {
        this.centerSortedFeatures = centerSortedFeatures;
    }

    public class FeatureRow {
        int start;
        int end;
        List<T> features;

        public FeatureRow() {
            this.features = new ArrayList(100);
        }

        public FeatureRow(List<T> features) {
            this.features = features;
        }

        public void addFeature(T feature) {
            if (features.isEmpty()) {
                this.start = getFeatureStartForPacking(feature);
            }
            features.add(feature);
            end = getFeatureEndForPacking(feature);
        }

        public List<T> getFeatures() {
            return features;
        }
    }
}
