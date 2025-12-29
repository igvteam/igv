package org.broad.igv.sam;

import org.broad.igv.feature.Range;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Convenience class for storing map from group -> list of alignments packed into rows
 * Can contain discontinuous regions
 * @author jacob
 * @date 2014-Jan-10
 */
public class PackedAlignments extends LinkedHashMap<String, List<Row>> {


    /**
     *  Set of ranges covered by this instance.
     */
    private List<? extends Range> ranges;

    PackedAlignments(List<? extends Range> ranges, Map<String, List<Row>> packedAlignments){
        super(packedAlignments);
        this.ranges = ranges;
    }

    /**
     * Calculate the total number of levels across all groups.
     * Since groups are stacked vertically, this is used to calculate height
     * @return
     */
    public int getNLevels() {
        int intervalNLevels = 0;
        for (List<Row> rows : this.values()) {
            intervalNLevels += rows.size();
        }
        return intervalNLevels;
    }

    /**
     * Whether any of the ranges contain this interval
     * @param chr
     * @param start
     * @param end
     * @return
     */
    public boolean contains(String chr, int start, int end) {
        for(Range range: this.ranges){
            if(range.contains(chr, start, end)){
                return true;
            }
        }
        return false;
    }

}
