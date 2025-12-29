package org.igv.sam;

import org.igv.event.DataLoadedEvent;
import org.igv.feature.genome.ChromAliasManager;
import org.igv.event.IGVEvent;
import org.igv.logging.*;
import org.igv.Globals;
import org.igv.event.IGVEventBus;
import org.igv.event.IGVEventObserver;
import org.igv.event.RefreshEvent;
import org.igv.feature.Range;
import org.igv.feature.genome.Genome;
import org.igv.prefs.IGVPreferences;
import org.igv.prefs.PreferencesManager;
import org.igv.sam.mods.BaseModificationKey;
import org.igv.sam.mods.BaseModificationSet;
import org.igv.sam.reader.AlignmentReader;
import org.igv.sam.reader.AlignmentReaderFactory;
import org.igv.track.Track;
import org.igv.ui.panel.FrameManager;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.util.ResourceLocator;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

import static org.igv.prefs.Constants.*;

/**
 * Manages data loading for a single alignment file.  Shared between alignment, coverage, and junction
 * tracks.
 */

public class AlignmentDataManager implements IGVEventObserver {

    private static Logger log = LogManager.getLogger(AlignmentDataManager.class);

    private final Genome genome;
    private ChromAliasManager chromAliasManager;
    private final AlignmentReader reader;
    private AlignmentTrack alignmentTrack;
    private CoverageTrack coverageTrack;
    private Set<Track> subscribedTracks;
    private List<AlignmentInterval> intervalCache;
    private ResourceLocator locator;

    private Set<String> sequenceNames;
    private AlignmentTileLoader loader;
    private Map<String, PEStats> peStats;
    private Set<BaseModificationKey> allBaseModificationKeys = new HashSet<>();
    private Set<String> simplexBaseModfications = new HashSet<>();

    private Range currentlyLoading;
    private AlignmentTrack.ExperimentType inferredType;

    Map<String, String> chrAliasCache = new HashMap<>();

    public AlignmentDataManager(ResourceLocator locator, Genome genome) throws IOException {
        this.locator = locator;
        // The time-gated limit for an AWS signed URL has expired, we need to re-sign the URL with the newly acquired
        // access token, otherwise we will face an Access Denied error. CheckReader() provides a very low overhead
        // mechanism to refresh expired presigned URLs.
        reader = AlignmentReaderFactory.getReader(locator);
        loader = new AlignmentTileLoader(reader);
        this.inferType();
        peStats = new HashMap();
        this.genome = genome;
        intervalCache = Collections.synchronizedList(new ArrayList<>());
        subscribedTracks = Collections.synchronizedSet(new HashSet<>());
        sequenceNames = new HashSet(getLoader().getSequenceNames());
        chromAliasManager = new ChromAliasManager(sequenceNames, genome);

        IGVEventBus.getInstance().subscribe(FrameManager.ChangeEvent.class, this);
        IGVEventBus.getInstance().subscribe(RefreshEvent.class, this);
    }

    public void receiveEvent(IGVEvent event) {
        if (event instanceof FrameManager.ChangeEvent) {
            trimCache();
        } else if (event instanceof RefreshEvent) {
            clear();
        } else {
            log.warn("Unknown event type: " + event.getClass());
        }
    }

    public void subscribe(Track track) {
        subscribedTracks.add(track);
    }

    public void unsubscribe(Track track) {
        subscribedTracks.remove(track);
        if (subscribedTracks.isEmpty()) {
            dispose();
            dumpAlignments();
            IGVEventBus.getInstance().unsubscribe(this);
        }
    }

    public AlignmentTileLoader getLoader() {
        return loader;
    }

    public ResourceLocator getLocator() {
        return locator;
    }

    public Map<String, PEStats> getPEStats() {
        return peStats;
    }

    /**
     * Return all base modfications seen in loaded alignments
     */
    public Set<BaseModificationKey> getAllBaseModificationKeys() {
        return allBaseModificationKeys;
    }

    public boolean isPairedEnd() {
        return getLoader().isPairedEnd();
    }

    public boolean hasYCTags() {
        return getLoader().hasYCTags();
    }

    public boolean hasIndex() {
        return getLoader().hasIndex();
    }

    public void setAlignmentTrack(AlignmentTrack alignmentTrack) {
        this.alignmentTrack = alignmentTrack;
    }

    public void setCoverageTrack(CoverageTrack coverageTrack) {
        this.coverageTrack = coverageTrack;
    }

    public CoverageTrack getCoverageTrack() {
        return coverageTrack;
    }

    public double getVisibilityWindow() {
        return getPreferences().getAsFloat(SAM_MAX_VISIBLE_RANGE) * 1000;
    }


    private IGVPreferences getPreferences() {
        String category = NULL_CATEGORY;
        AlignmentTrack.ExperimentType experimentType = getExperimentType();
        if (experimentType == AlignmentTrack.ExperimentType.RNA) {
            category = RNA;
        } else if (experimentType == AlignmentTrack.ExperimentType.THIRD_GEN) {
            category = THIRD_GEN;
        }
        return PreferencesManager.getPreferences(category);
    }


    /**
     * The set of sequences found in the file.
     * May be null
     *
     * @return
     */
    public List<String> getSequenceNames() throws IOException {
        return getLoader().getSequenceNames();
    }

    public AlignmentInterval getLoadedInterval(ReferenceFrame frame) {
        return getLoadedInterval(frame, false);
    }

    public AlignmentInterval getLoadedInterval(ReferenceFrame frame, boolean includeOverlaps) {
        // Search for interval completely containing reference frame region
        for (AlignmentInterval interval : intervalCache) {
            if (interval.contains(frame.getCurrentRange())) {
                return interval;
            }
        }
        // No interval contains entire region of frame, look for intervael with at least some overlap
        if (includeOverlaps) {
            for (AlignmentInterval interval : intervalCache) {
                if (interval.overlaps(frame.getCurrentRange())) {
                    return interval;
                }
            }
        }
        return null;
    }

    public void setViewAsPairs(boolean option, AlignmentTrack.RenderOptions renderOptions, Track.DisplayMode displayMode) {
        if (option == renderOptions.isViewPairs()) {
            return;
        }
        renderOptions.setViewPairs(option);
        packAlignments(renderOptions, displayMode);
    }

    /**
     * Repack currently loaded alignments across frames
     * All relevant intervals must be loaded
     *
     * @param renderOptions
     * @return Whether repacking was performed
     */
    void packAlignments(AlignmentTrack.RenderOptions renderOptions, Track.DisplayMode displayMode) {
        for (AlignmentInterval interval : intervalCache) {
            interval.packAlignments(renderOptions, displayMode);
        }
    }


    public boolean isLoaded(ReferenceFrame frame) {
        return getLoadedInterval(frame) != null;
    }

    public void load(ReferenceFrame frame,
                     AlignmentTrack.RenderOptions renderOptions,
                     Track.DisplayMode displayMode,
                     boolean expandEnds) {

        if (frame.getChrName().equals(Globals.CHR_ALL) || (frame.getEnd() - frame.getOrigin()) > getVisibilityWindow())
            return; // should not happen

        if (isLoaded(frame)) {
            return;  // Already loaded
        }

        Range range = frame.getCurrentRange();
        final String chr = frame.getChrName();
        final int start = range.getStart();
        final int end = range.getEnd();
        int adjustedStart = start;
        int adjustedEnd = end;
        Range adjustedRange = new Range(chr, start, end);

        if (currentlyLoading != null && currentlyLoading.contains(adjustedRange)) {
            return;  // Already loading
        }
        try {
            currentlyLoading = adjustedRange;

            // Expand the interval by the lesser of  +/- a 2 screens, or max visible range
            int windowSize = Math.min(4 * (end - start), PreferencesManager.getPreferences().getAsInt(SAM_MAX_VISIBLE_RANGE) * 1000);
            int center = start + (end - start) / 2;     // Be careful how you calculate this -- potential overflow for large chromosomes
            int expand = Math.min(Integer.MAX_VALUE - center, Math.max(end - start, windowSize / 2));

            if (expandEnds) {
                adjustedStart = Math.max(0, Math.min(start, center - expand));
                adjustedEnd = Math.max(end, center + expand);
                if (adjustedEnd < 0) {
                    adjustedEnd = Integer.MAX_VALUE;  // Overflow
                }
            }

            AlignmentInterval loadedInterval = loadInterval(chr, adjustedStart, adjustedEnd, renderOptions, frame);

            trimCache();

            intervalCache.add(loadedInterval);

            loadedInterval.packAlignments(renderOptions, displayMode);

            IGVEventBus.getInstance().post(new DataLoadedEvent(frame));

        } finally {
            currentlyLoading = null;
        }
    }

    /**
     * Remove out-of-view intervals from the cache.  This is O(N) where N = #frames X #intervals.   It is assumed
     * that N is small
     */
    private void trimCache() {
        List<AlignmentInterval> trimmedIntervals = new ArrayList<>();
        for (AlignmentInterval interval : intervalCache) {
            if (intervalInView(interval)) {
                trimmedIntervals.add(interval);
            }
        }
        intervalCache = trimmedIntervals;
    }


    private boolean intervalInView(AlignmentInterval interval) {

        for (ReferenceFrame frame : FrameManager.getFrames()) {
            if (interval.contains(frame.getCurrentRange())) {
                return true;
            }
        }
        return false;
    }


    AlignmentInterval loadInterval(String chr, int start, int end, AlignmentTrack.RenderOptions renderOptions, ReferenceFrame frame) {

        final String seqName;
        if (sequenceNames.contains(chr)) {
            seqName = chr;
        } else if (chrAliasCache.containsKey(chr)) {
            seqName = chrAliasCache.get(chr);
        } else {
            seqName = chromAliasManager.getAliasName(chr);
            chrAliasCache.put(chr, seqName);   // This might be null, but record to prevent further attempts
        }

        if (seqName == null) {
            // No alignments with this chr name -- return empty interval
            return new AlignmentInterval(chr, start, end);
        }

        DownsampleOptions downsampleOptions = new DownsampleOptions();
        final AlignmentTrack.BisulfiteContext bisulfiteContext =
                renderOptions != null ? renderOptions.bisulfiteContext : null;
        SpliceJunctionHelper spliceJunctionHelper = new SpliceJunctionHelper();

        AlignmentTileLoader.AlignmentTile t = getLoader().loadTile(seqName, start, end, spliceJunctionHelper,
                downsampleOptions, peStats, bisulfiteContext, renderOptions);
        List<Alignment> alignments = t.getAlignments();
        List<DownsampledInterval> downsampledIntervals = t.getDownsampledIntervals();
        this.updateBaseModfications(alignments);
        return new AlignmentInterval(chr, start, end, alignments, t.getCounts(), spliceJunctionHelper, downsampledIntervals, frame);
    }

    private void updateBaseModfications(List<Alignment> alignments) {

        for (Alignment a : alignments) {
            List<BaseModificationSet> bmSets = a.getBaseModificationSets();
            if (bmSets != null) {
                for (BaseModificationSet bms : bmSets) {
                    allBaseModificationKeys.add(BaseModificationKey.getKey(bms.getBase(), bms.getStrand(), bms.getModification()));
                }
            }
        }

        // Search for simplex modifications (single strand read, e.g. C+m with no G-m)
        Set<String> minusStranMods = allBaseModificationKeys.stream()
                .filter(key -> key.getStrand() == '-')
                .map(key -> key.getModification())
                .collect(Collectors.toSet());
        for (BaseModificationKey key : allBaseModificationKeys) {
            if (key.getStrand() == '+' && !minusStranMods.contains(key.getModification()))
                simplexBaseModfications.add(key.getModification());
            simplexBaseModfications.add("NONE_" + key.getCanonicalBase());  // Mix of simplex & duplex keys for same base not supported.
        }
    }

    public AlignmentTrack.ExperimentType inferType() {

        if (this.inferredType == null) {
            try {
                ReadStats readStats = new ReadStats();
                List<Alignment> sample = AlignmentUtils.firstAlignments(reader, 100);
                for (Alignment a : sample) {
                    readStats.addAlignment(a);
                }
                inferredType = readStats.inferType();

                if (inferredType == AlignmentTrack.ExperimentType.THIRD_GEN) {
                    return inferredType;
                } else {
                    // Get a larger sample to distinguish RNA-Seq
                    readStats = new ReadStats();
                    sample = AlignmentUtils.firstAlignments(reader, 2000);
                    for (Alignment a : sample) {
                        readStats.addAlignment(a);
                    }
                    inferredType = readStats.inferType();
                    return inferredType;
                }
            } catch (Exception e) {
                // Ignore errors, inferring alignment type is not essential.
                log.error("Error inferring alignment type", e);
                inferredType = AlignmentTrack.ExperimentType.UNKOWN;
            }
        }

        return inferredType;
    }

    private AlignmentTrack.ExperimentType getExperimentType() {
        return alignmentTrack == null ? null : alignmentTrack.getExperimentType();
    }


    public PackedAlignments getGroups(AlignmentInterval interval, AlignmentTrack.RenderOptions renderOptions) {
        //AlignmentInterval interval = getLoadedInterval(context.getReferenceFrame());
        if (interval != null) {
            return interval.getPackedAlignments();
        } else {
            return null;
        }
    }

    public void clear() {
        intervalCache.clear();
    }

    public void dumpAlignments() {
        for (AlignmentInterval interval : intervalCache) {
            interval.dumpAlignments();
        }
    }

    /**
     * Find the first loaded interval for the specified chromosome and genomic {@code positon},
     * return the grouped alignments
     *
     * @param position
     * @param frame
     * @return alignmentRows, grouped and ordered by key
     */
    public PackedAlignments getGroupedAlignmentsContaining(double position, ReferenceFrame frame) {
        String chr = frame.getChrName();
        int start = (int) position;
        int end = start + 1;

        AlignmentInterval interval = getLoadedInterval(frame);
        if (interval == null) {
            return null;
        } else {
            PackedAlignments packedAlignments = interval.getPackedAlignments();
            if (packedAlignments != null && packedAlignments.contains(chr, start, end)) {
                return packedAlignments;
            } else {
                return null;
            }
        }
    }

    public int getNLevels() {
        int nLevels = 0;

        for (AlignmentInterval interval : intervalCache) {
            PackedAlignments packedAlignments = interval.getPackedAlignments();
            if (packedAlignments != null) {
                int intervalNLevels = packedAlignments.getNLevels();
                nLevels = Math.max(nLevels, intervalNLevels);
            }
        }
        return nLevels;
    }

    /**
     * Get the maximum group count among all the loaded intervals.  Normally there is one interval, but there
     * can be multiple if viewing split screen.
     */
    public int getMaxGroupCount() {
        int groupCount = 0;

        for (AlignmentInterval interval : intervalCache) {
            if (interval != null) {  // Not sure how this happens but it does
                PackedAlignments packedAlignments = interval.getPackedAlignments();
                if (packedAlignments != null) {
                    groupCount = Math.max(groupCount, packedAlignments.size());
                }
            }
        }
        return groupCount;
    }


    private void dispose() {
        if (loader != null) {
            try {
                loader.close();
            } catch (IOException ex) {
                log.error("Error closing AlignmentQueryReader. ", ex);
            }
        }
    }

    public void updatePEStats(AlignmentTrack.RenderOptions renderOptions) {
        if (this.peStats != null) {
            for (PEStats stats : peStats.values()) {
                stats.computeInsertSize(renderOptions.getMinInsertSizePercentile(), renderOptions.getMaxInsertSizePercentile());
            }
        }
    }

    public void alleleThresholdChanged() {
        coverageTrack.setSnpThreshold(PreferencesManager.getPreferences().getAsFloat(SAM_ALLELE_THRESHOLD));
    }

    public boolean isPhased() {
        return getLoader().isPhased();
    }

    public Collection<AlignmentInterval> getLoadedIntervals() {
        return intervalCache;
    }

    public Set<String> getSimplexBaseModifications() {
        return simplexBaseModfications;
    }

    public static class DownsampleOptions {
        private boolean downsample;
        private int sampleWindowSize;
        private int maxReadCount;

        public DownsampleOptions() {
            IGVPreferences prefs = PreferencesManager.getPreferences();
            init(prefs.getAsBoolean(SAM_DOWNSAMPLE_READS),
                    prefs.getAsInt(SAM_SAMPLING_WINDOW),
                    prefs.getAsInt(SAM_SAMPLING_COUNT));
        }

        DownsampleOptions(boolean downsample, int sampleWindowSize, int maxReadCount) {
            init(downsample, sampleWindowSize, maxReadCount);
        }

        private void init(boolean downsample, int sampleWindowSize, int maxReadCount) {
            this.downsample = downsample;
            this.sampleWindowSize = sampleWindowSize;
            this.maxReadCount = maxReadCount;
        }

        public boolean isDownsample() {
            return downsample;
        }

        public int getSampleWindowSize() {
            return sampleWindowSize;
        }

        public int getMaxReadCount() {
            return maxReadCount;
        }

    }

}

