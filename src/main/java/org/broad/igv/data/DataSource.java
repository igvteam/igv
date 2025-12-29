/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.data;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;

import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * @author jrobinso
 */
public interface DataSource {

    double getDataMax();

    double getDataMin();

    List<LocusScore> getSummaryScoresForRange(
            String chr,
            int startLocation,
            int endLocation,
            int zoom);

    TrackType getTrackType();

    void setWindowFunction(WindowFunction statType);

    boolean isLogNormalized();

    public WindowFunction getWindowFunction();

    public default Collection<WindowFunction> getAvailableWindowFunctions() {
        return Collections.EMPTY_LIST;
    }

    default boolean isIndexable() {
        return true;
    }

    /**
     * Release any resources (file handles, etc)
     */
    void dispose();
}
