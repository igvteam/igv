/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

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
