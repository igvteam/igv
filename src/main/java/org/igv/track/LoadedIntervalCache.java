package org.igv.track;

import org.igv.event.IGVEvent;
import org.igv.event.IGVEventBus;
import org.igv.event.IGVEventObserver;
import org.igv.feature.Chromosome;
import org.igv.feature.Locus;
import org.igv.feature.genome.Genome;
import org.igv.feature.genome.GenomeManager;
import org.igv.ui.panel.FrameManager;
import org.igv.ui.panel.ReferenceFrame;

import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * Per-frame cache of loaded data intervals, keyed by frame name, for tracks that load one interval per reference
 * frame.  Entries for frames that no longer exist are purged on FrameManager.ChangeEvent.
 */
public class LoadedIntervalCache<T> implements IGVEventObserver {

    private Map<String, LoadedDataInterval<T>> cache = Collections.synchronizedMap(new HashMap<>());

    public LoadedIntervalCache() {
        IGVEventBus.getInstance().subscribe(FrameManager.ChangeEvent.class, this);
    }

    public LoadedDataInterval<T> get(String frameName) {
        return cache.get(frameName);
    }

    public void put(String frameName, LoadedDataInterval<T> interval) {
        cache.put(frameName, interval);
    }

    public Collection<LoadedDataInterval<T>> values() {
        return cache.values();
    }

    public void clear() {
        cache.clear();
    }

    /**
     * True if the interval loaded for the frame covers the frame's current chromosome, range and zoom.
     */
    public boolean contains(ReferenceFrame frame) {
        LoadedDataInterval<T> interval = cache.get(frame.getName());
        return interval != null && interval.contains(frame);
    }

    public void unsubscribe() {
        IGVEventBus.getInstance().unsubscribe(this);
    }

    @Override
    public void receiveEvent(IGVEvent event) {
        if (event instanceof FrameManager.ChangeEvent) {
            // Purge any cached data not associated with current frames
            Collection<ReferenceFrame> frames = ((FrameManager.ChangeEvent) event).frames();
            Map<String, LoadedDataInterval<T>> newCache = Collections.synchronizedMap(new HashMap<>());
            for (ReferenceFrame f : frames) {
                newCache.put(f.getName(), cache.get(f.getName()));
            }
            cache = newCache;
        }
    }

    /**
     * The interval to load for a frame: the view expanded by 50% on each side (for panning) and clamped to the
     * chromosome, unless in a multi-locus mode with "lots" of frames.  The chromosome name is the frame's.
     */
    public static Locus expandedInterval(ReferenceFrame frame) {

        String chr = frame.getChrName();
        int start = (int) frame.getOrigin();
        int end = (int) frame.getEnd() + 1;
        int maxEnd = end;
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        if (genome != null) {
            Chromosome c = genome.getChromosome(chr);
            if (c != null) maxEnd = Math.max(c.getLength(), end);
        }

        boolean multiLocus = (FrameManager.getFrames().size() > 4);
        int delta = multiLocus ? 1 : (end - start) / 2;
        int expandedStart = Math.max(0, start - delta);
        int expandedEnd = Math.min(maxEnd, end + delta);
        if (expandedEnd < 0) {
            // overflow
            expandedEnd = Integer.MAX_VALUE;
        }
        return new Locus(chr, expandedStart, expandedEnd);
    }
}
