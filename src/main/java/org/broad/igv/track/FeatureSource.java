package org.broad.igv.track;

import htsjdk.tribble.Feature;
import htsjdk.tribble.NamedFeature;
import org.broad.igv.bedpe.BedPE;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * A queryable source for Features.  Used by FeatureTrack
 *
 * author: jrobinso
 * Date: Jan 31, 2010
 */
public interface FeatureSource<T extends Feature> {


    /**
     * Return an iterator over all features that overlap the interval.  The coordinates are in the "UCSC" convention,
     * that is first base is zero and the interval is not inclusive of the end position,  i.e. an interval
     * spanning the first base is 0 - 1.
     *
     * @param chr
     * @param start
     * @param end
     * @return
     */
    Iterator<T> getFeatures(String chr, int start, int end) throws IOException;


    /**
     * Return a list of coverage values spanning the given interval.  This can be null if coverage is not known
     * or relevant.
     *
     * @param chr
     * @param start
     * @param end
     * @param zoom   the zoom level
     * @return
     */
    default List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {
        return null;
    }

    /**
     * The featureWindowSize is the genomic interval, in bases, at which the associated feature track should start
     * rendering features.   When zoomed out beyond this interval coverage is shown, if available, or a message to
     * zoom in.
     *
     * TODO:  It seems slightly odd to have this property here, but on the other hand FeatureSource objects have
     * TODO:  access to the information required to estimate this value.  Still, consider removing this from the interface.
     *
     * @return the feature window threshold in base-pairs
     */
    default int getFeatureWindowSize() {
        return -1;
    }

    default void close() {
      // Do nothing by default
    };

    default Object getHeader() throws IOException {
        return null;
    }

    /**
     * Return true if the source can be searched for a feature by name
     *
     * @return
     */
    default boolean isSearchable() {
        return false;
    }

    default List<NamedFeature> search(String name) {
        return null;
    }
}
