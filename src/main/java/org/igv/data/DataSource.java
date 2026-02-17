
package org.igv.data;

import org.igv.feature.LocusScore;
import org.igv.track.DataType;
import org.igv.track.WindowFunction;

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

    DataType getDataType();

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
