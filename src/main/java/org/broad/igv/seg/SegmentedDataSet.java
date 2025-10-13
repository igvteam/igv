
package org.broad.igv.seg;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.TrackType;

import java.util.List;

/**
 * @author jrobinso
 */
public interface SegmentedDataSet {

    double getDataMax(String chr);

    double getDataMin(String chr);

    List<String> getSampleNames();

    List<LocusScore> getSegments(String heading, String chr);

    TrackType getType();

    List<LocusScore> getWholeGenomeScores(String heading);

    boolean isLogNormalized();

}
