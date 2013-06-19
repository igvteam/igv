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

/**
 * @author jacob
 * @date 2013-Jun-17
 */
public interface IAlignmentDataManager {

    /**
     * Return the loaded interval for the specified frame (by name).  Note this can be null if the interval isn't loaded
     * yet.
     *
     * @param frameName
     * @return
     */
    public AlignmentInterval getLoadedInterval(String frameName);

    SpliceJunctionHelper.LoadOptions getSpliceJunctionLoadOptions();

    Collection<AlignmentInterval> getAllLoadedIntervals();

    void setMinJunctionCoverage(int newMinJunctionCoverage);
}
