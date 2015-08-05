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

import java.util.Collection;

/**
 * Implementation of IAlignmentDataManager which simply caches loaded intervals,
 * and returns them on request
 *
 * @author jacob
 * @date 2013-Jun-17
 */
public class MemoryAlignmentDataManager implements IAlignmentDataManager {


    private SpliceJunctionHelper.LoadOptions loadOptions;

    private PositionCache<AlignmentInterval> loadedIntervalCache = new PositionCache<AlignmentInterval>();

    public MemoryAlignmentDataManager(AlignmentDataManager alignmentDataManager, SpliceJunctionHelper.LoadOptions loadOptions) {
        this.loadOptions = loadOptions;
        this.loadedIntervalCache = new PositionCache(alignmentDataManager.getCache());
    }

    @Override
    public AlignmentInterval getLoadedInterval(Range range){
        return this.loadedIntervalCache.getForRange(range);
    }

    @Override
    public SpliceJunctionHelper.LoadOptions getSpliceJunctionLoadOptions() {
        return loadOptions;
    }

    @Override
    public void setMinJunctionCoverage(int minJunctionCoverage) {
        this.loadOptions = new SpliceJunctionHelper.LoadOptions(minJunctionCoverage, this.loadOptions.minReadFlankingWidth);
        for (AlignmentInterval interval : getLoadedIntervals()) {
            interval.getSpliceJunctionHelper().setLoadOptions(this.loadOptions);
        }
    }

    @Override
    public Collection<AlignmentInterval> getLoadedIntervals() {
        return loadedIntervalCache.values();
    }


}
