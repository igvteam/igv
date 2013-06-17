/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.sam;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

/**
 * Implementation of IAlignmentDataManager which simply caches loaded intervals,
 * and returns them on request
 *
 * @author jacob
 * @date 2013-Jun-17
 */
public class MemoryAlignmentDataManager implements IAlignmentDataManager {


    private SpliceJunctionHelper.LoadOptions loadOptions;

    private Map<String, AlignmentInterval> loadedIntervalMap = new HashMap<String, AlignmentInterval>();

    public MemoryAlignmentDataManager(AlignmentDataManager alignmentDataManager, SpliceJunctionHelper.LoadOptions loadOptions) {
        this.loadOptions = loadOptions;
        for(String intervalName: alignmentDataManager.getLoadedIntervalNames()){
            this.loadedIntervalMap.put(intervalName, alignmentDataManager.getLoadedInterval(intervalName));
        }

    }

    @Override
    public AlignmentInterval getLoadedInterval(String intervalName) {
        return this.loadedIntervalMap.get(intervalName);
    }

    @Override
    public SpliceJunctionHelper.LoadOptions getSpliceJunctionLoadOptions() {
        return loadOptions;
    }

    @Override
    public void setMinJunctionCoverage(int minJunctionCoverage) {
        this.loadOptions = new SpliceJunctionHelper.LoadOptions(minJunctionCoverage, this.loadOptions.minReadFlankingWidth);
        for (AlignmentInterval interval : getAllLoadedIntervals()) {
            interval.getSpliceJunctionHelper().setLoadOptions(this.loadOptions);
        }
    }

    @Override
    public Collection<AlignmentInterval> getAllLoadedIntervals() {
        return loadedIntervalMap.values();
    }


}
