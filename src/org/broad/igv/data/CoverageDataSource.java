/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.data;

/**
 * @author Fabien Campagne
 *         Date: 6/10/11
 *         Time: 4:45 PM
 */
public interface CoverageDataSource extends DataSource {
    /**
     * The filename that contains the coverage data. Used to persist the state of coverage tracks.
     * @return a filename.
     */
    String getPath();

    /**
     * Tell the coverage source to normalize coverage by some appropriate normalization method.
     * @param normalize True if normalization should be performed, false otherwise.
     */
    void setNormalize(boolean normalize);

    /**
     * Whether the source is performing
     * some normalization operation
     */
    boolean getNormalize();
}
