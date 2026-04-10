package org.igv.sample;

import org.igv.track.RegionScoreType;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public interface SamplesTrack {

    /**
     * Return the number of samples in this track.  For most tracks this will be 1.
     *
     * @return
     */
    default int sampleCount() {
        return 1;
    }

    default List<SampleGroup> getSampleGroups() {
        return Collections.EMPTY_LIST;
    }

    default int getSampleHeight() {
        return 15;
    }

    default int getSampleOffset() {
        return 0;
    }

    default List<String> getSampleNames() {
        return Collections.emptyList();
    }

    default void sortSamples(Comparator<String> comparator) {
    }

    default void sortSamplesByValue(String chr, int start, int end, RegionScoreType type) {
        // no op, override in subclass if needed
    }

    default void filterSamples(SampleFilter sampleFilter) {
        // no op, override in subclass if needed
    }

}
