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

import org.broad.igv.feature.Range;

import java.util.Collection;

/**
 * @author jacob
 * @date 2013-Jun-17
 */
public interface IAlignmentDataManager {

    /**
     * Return the loaded interval for the specified range.
     * This can be null if it's not loaded (or if loading is in progress)
     * @param range
     * @return
     */
    public AlignmentInterval getLoadedInterval(Range range);

    SpliceJunctionHelper.LoadOptions getSpliceJunctionLoadOptions();

    Collection<AlignmentInterval> getLoadedIntervals();

    void setMinJunctionCoverage(int newMinJunctionCoverage);
}
