package org.igv.feature;

import org.igv.feature.genome.ChromAlias;

import java.awt.*;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Interface for features in IGV annotation tracks  (FeatureTrack and derived classes).
 */

public interface IGVFeature extends LocusScore, IGVNamedFeature {

    default String getIdentifier() {
        return null;
    }

    default Strand getStrand() {
        return Strand.NONE;
    }

    default int getLength() {
        return getEnd() - getStart();
    }

    default List<String> getAttributeKeys() {
        return Collections.EMPTY_LIST;
    }

    default String getAttribute(String key) {
        return null;
    }

    default void removeAttribute(String key) {}

    /**
     * Return true if the given feature is completely contained within the bounds of this
     * feature. amd is on the same strand.
     * <p/>
     *
     * @param feature
     * @return
     */
    default boolean contains(IGVFeature feature) {
        if (feature == null) {
            return false;
        }
        if (!this.getChr().equals(feature.getChr()) ||
                this.getStrand() != feature.getStrand()) {
            return false;
        }
        if ((feature.getStart() >= this.getStart()) && (feature.getEnd() <= this.getEnd())) {
            return true;
        } else {
            return false;
        }
    }

    default boolean contains(double location) {
        return location >= getStart() && location < getEnd();
    }

    default List<Exon> getExons() {
        return null;
    }

    public Color getColor();

    default String getURL() {
        return null;
    }

   default Map<String, String> getAttributes() {
        return Collections.emptyMap();
    }
}
