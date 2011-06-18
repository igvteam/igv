/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.h5;

import java.util.List;

/**
 * An HDF55 reader handles read access to an HDF5 file,  either local or
 * remote.
 *
 * @author jrobinso
 */
public interface HDF5Reader {

    public void closeFile();

    /**
     * Open an HDF5 dataset in the given file.
     *
     * @param datasetName
     * @return a handle to the dataset
     */
    int openDataset(String datasetName);

    /**
     * @param datasetId
     */
    void closeDataset(int datasetId);

    /**
     * Open an HDF5 data group.
     *
     * @param groupPath
     * @return
     */
    int openGroup(String groupPath);

    /**
     * Close an hdf5 group
     *
     * @param groupId
     */
    public void closeGroup(int groupId);

    /**
     * Return value of the integer attribute for the specified group
     *
     * @param groupId
     * @param attrName
     * @return
     *
     */
    public int readIntegerAttribute(int groupId, String attrName) throws ObjectNotFoundException;

    /**
     * Return value of the string attribute for the specified group
     *
     * @param groupId
     * @param attrName
     * @return
     *
     */
    public String readStringAttribute(int groupId, String attrName);

    /**
     * Return value of the double attribute for the specified group
     *
     * @param groupId
     * @param attrName
     * @return
     */
    public double readDoubleAttribute(int groupId, String attrName);

    /**
     * Read the entire dataset as an array of floats.
     *
     * @param datasetId
     * @return
     */
    float[] readAllFloats(int datasetId);

    /**
     * Read the entire dataset as an array of ints.
     *
     * @param datasetId
     * @return
     */
    int[] readAllInts(int datasetId);

    /**
     * Read the entire dataset as an array of strings.
     *
     * @param datasetId
     * @return
     */
    List<String> readAllStrings(int datasetId);

    /**
     * Read a section of the dataset as an array of floats.
     *
     * @param datasetId
     * @param fromIndex
     * @param toIndex
     * @return
     */
    public float[] readFloats(int datasetId, int fromIndex, int toIndex);

    /**
     * Read a section of the dataset as an array of ints.
     *
     * @param datasetId
     * @param fromIndex
     * @param toIndex
     * @return
     */
    public int[] readInts(int datasetId, int fromIndex, int toIndex);

    /**
     * Read a section of the dataset as a 2-D array of floats.
     *
     * @param fromIndex Start column index
     * @param toIndex   End column index
     * @return
     */
    public float[][] readDataSlice(int datasetId, int fromIndex, int toIndex);
}
