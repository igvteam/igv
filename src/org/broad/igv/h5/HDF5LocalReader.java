/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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


/*
* IGVTestReadH5Lib.java
*
* Tests reading an IGV H5 file using the low level C wrapper library.
* This file maps java functions as directly as possible to the underlying
* C functions.  It is much faster than the object API.
*/
package org.broad.igv.h5;

//~--- non-JDK imports --------------------------------------------------------

import ncsa.hdf.hdf5lib.H5;
import ncsa.hdf.hdf5lib.HDF5Constants;
import ncsa.hdf.hdf5lib.exceptions.HDF5AttributeException;
import ncsa.hdf.hdf5lib.exceptions.HDF5Exception;
import ncsa.hdf.hdf5lib.exceptions.HDF5LibraryException;
import org.apache.log4j.Logger;
import org.broad.igv.util.ObjectCache;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 */
public class HDF5LocalReader implements HDF5Reader {

    private static Logger log = Logger.getLogger(HDF5LocalReader.class);

    /**
     * Field description
     */
    public static int RDONLY = ncsa.hdf.hdf5lib.HDF5Constants.H5F_ACC_RDONLY;

    static ObjectCache<Integer, HDF5EntityCache> entityCaches = new ObjectCache();

    private int fileId = -1;

    /**
     * Constructs ...
     *
     * @param path
     */
    public HDF5LocalReader(String path) {
        try {
            this.fileId = openFile(path, RDONLY);
        }
        catch (FileNotFoundException ex) {
            ex.printStackTrace();
        }
    }

    /**
     * Open an HDF 5 and return a hnadle to the file.  The handle
     * should be stored by the caller and used for further
     * file operaions.
     *
     * @param path a unique identifier for the file.  Will be the file path for a
     *             local reader, but could be something else for a remote reader.
     * @return handle to the file (basically a C pointer)
     * @throws java.io.FileNotFoundException
     */
    private int openFile(String path, int mode) throws FileNotFoundException {
        try {
            fileId = H5.H5Fopen(path, mode, ncsa.hdf.hdf5lib.HDF5Constants.H5P_DEFAULT);
            entityCaches.put(fileId, new HDF5EntityCache(fileId));
            return fileId;
        }
        catch (HDF5LibraryException ex) {
            String msg = "Error opening HDF5 file: " + path;
            log.error(msg, ex);
            throw new DataAccessException(msg, ex);
        }
    }

    /**
     * Close an hdf5 file
     */
    public void closeFile() {
        try {
            if (fileId > 0) {
                H5.H5Fclose(fileId);
                HDF5EntityCache cache = entityCaches.get(fileId);
                if (cache != null) {
                    cache.closeAllEntities();
                }
                fileId = -1;
            }
        }
        catch (HDF5LibraryException ex) {
            String msg = "Error closing HDF5 file";
            log.error(msg, ex);
            throw new DataAccessException(msg, ex);
        }
    }

    /**
     * Method description
     *
     * @param dsName
     * @return
     */
    public int openDataset(String dsName) {
        HDF5EntityCache cache = entityCaches.get(fileId);
        if (cache == null) {
            cache = new HDF5EntityCache(fileId);
            entityCaches.put(fileId, cache);
        }
        return cache.openDataset(dsName);
    }

    /**
     * Close an hdf5 dataset
     *
     * @param datasetId
     */
    public void closeDataset(int datasetId) {

        // Closing now managed by cache
    }

    /**
     * Method description
     *
     * @param name
     * @return
     */
    public int openGroup(String name) {
        try {
            return H5.H5Gopen(fileId, name);
        }
        catch (ncsa.hdf.hdf5lib.exceptions.HDF5SymbolTableException ex) {
            String msg = "Error opening group: " + name;
            if (log.isDebugEnabled()) {
                log.debug(msg, ex);
            }
            throw new ObjectNotFoundException(msg);
        }
        catch (HDF5LibraryException ex) {
            String msg = "Error opening group: " + name;
            log.error(msg, ex);
            throw new DataAccessException(msg, ex);
        }

    }

    /**
     * Close an hdf5 group
     *
     * @param groupId
     */
    public void closeGroup(int groupId) {
        try {
            H5.H5Gclose(groupId);
        }
        catch (HDF5LibraryException ex) {
            String msg = "Error closing HDF5 group";
            log.error(msg, ex);
            throw new DataAccessException(msg, ex);
        }
    }

    /**
     * Returns the HDF typeId and array length for an array.  Info is returned as
     * a 2 element int array.  Yes, a class could be used but seems excessive
     * for this private method.
     */
    private int getTypeForArray(Object array) {

        if ((array instanceof Integer) || (array instanceof int[]) || (array instanceof int[][])) {
            return HDF5Constants.H5T_NATIVE_INT;
        }

        if ((array instanceof Long) || (array instanceof long[]) || (array instanceof long[][])) {
            return HDF5Constants.H5T_NATIVE_LLONG;
        }

        if ((array instanceof Float) || (array instanceof float[]) || (array instanceof float[][])) {
            return HDF5Constants.H5T_NATIVE_FLOAT;
        }

        if ((array instanceof Number) || (array instanceof double[])
                || (array instanceof double[][])) {
            return HDF5Constants.H5T_NATIVE_DOUBLE;
        }

        if (array instanceof char[]) {
            return createStringDatatype(1);
        }

        String msg = "Error: No HDF Type for: " + array.getClass().getName();
        log.error(msg);
        throw new DataAccessException("No HDF Type for: " + array.getClass().getName());

    }

    /**
     * Method description
     *
     * @param strLength
     * @return
     */
    public int createStringDatatype(int strLength) {
        try {
            int tid = H5.H5Tcopy(HDF5Constants.H5T_C_S1);
            H5.H5Tset_size(tid, strLength);
            return tid;
        }
        catch (HDF5LibraryException ex) {
            String msg = "Error creating string datatype. Length = " + strLength;
            log.error(msg, ex);
            throw new DataAccessException(msg, ex);
        }

    }

    /**
     * Read an interval from a single data row from the 2-D array.
     *
     * @param datasetId
     * @param rowId     The row number (first row is zero)
     * @param fromIndex Start column index
     * @param toIndex   End column index
     * @return
     */
    public float[] readDataRowSlice(int datasetId, int rowId, int fromIndex, int toIndex) {
        try {
            int dataspace = getDataspace(datasetId);
            int nPts = toIndex - fromIndex + 1;

            long[] offset = new long[2];
            long[] count = new long[2];
            offset[0] = rowId;
            offset[1] = fromIndex;
            count[0] = 1;
            count[1] = nPts;
            int status1 = H5.H5Sselect_hyperslab(dataspace, HDF5Constants.H5S_SELECT_SET, offset,
                    null, count, null);

            // Step 2.   Create a memory dataspace
            long[] memDims = new long[1];
            memDims[0] = nPts;
            int memspace = H5.H5Screate_simple(1, memDims, null);

            // Step 4.  Read the data
            float[] partialData = new float[nPts];
            H5.H5Dread(datasetId, HDF5Constants.H5T_NATIVE_FLOAT, memspace, dataspace,
                    HDF5Constants.H5P_DEFAULT, partialData);

            return partialData;

        }
        catch (HDF5Exception ex) {
            String msg = "Error reading data row";
            log.error(msg);
            throw new DataAccessException(msg);
        }

    }

    /**
     * Read a slice defined by fromIndex and toIndex from a 2-D array
     *
     * @param datasetId
     * @param fromIndex The start column index
     * @param toIndex   The end column index
     * @return
     */
    public float[][] readDataSlice(int datasetId, int fromIndex, int toIndex) {
        try {
            int dataspace = getDataspace(datasetId);
            int nPts = toIndex - fromIndex + 1;

            long[] dims = getDims(datasetId);
            int nRows = (int) dims[0];

            long[] offset = new long[2];
            long[] count = new long[2];
            offset[0] = 0;
            offset[1] = fromIndex;
            count[0] = nRows;
            count[1] = nPts;
            int status1 = H5.H5Sselect_hyperslab(dataspace, HDF5Constants.H5S_SELECT_SET, offset,
                    null, count, null);

            // Step 2.   Create a memory dataspace
            long[] memDims = new long[2];
            memDims[0] = nRows;
            memDims[1] = nPts;
            int memspace = H5.H5Screate_simple(2, memDims, null);

            // Step 4.  Read the data
            float[][] partialData = new float[nRows][nPts];
            H5.H5Dread(datasetId, HDF5Constants.H5T_NATIVE_FLOAT, memspace, dataspace,
                    HDF5Constants.H5P_DEFAULT, partialData);

            return partialData;

        }
        catch (HDF5Exception ex) {
            String msg = "Error reading data";
            log.error(msg, ex);
            throw new DataAccessException(msg, ex);
        }

    }

    /**
     * @param groupId
     * @param name
     * @return
     */
    public double readDoubleAttribute(int groupId, String name) {

        try {
            int attrId = H5.H5Aopen_name(groupId, name);
            double[] buffer = new double[1];
            H5.H5Aread(attrId, HDF5Constants.H5T_NATIVE_DOUBLE, buffer);
            H5.H5Aclose(attrId);
            return buffer[0];
        }
        catch (HDF5AttributeException e) {
            throw new ObjectNotFoundException("Attribute not found: " + name);
        }
        catch (HDF5Exception ex) {
            String msg = "Error reading attribute from dataset";
            log.error(msg, ex);
            throw new DataAccessException(msg, ex);
        }

    }

    /**
     * @param groupId
     * @param name
     * @return
     * @throws ObjectNotFoundException
     */
    public int readIntegerAttribute(int groupId, String name) throws ObjectNotFoundException {


        try {
            int attrId = H5.H5Aopen_name(groupId, name);
            int[] buffer = new int[1];
            H5.H5Aread(attrId, HDF5Constants.H5T_NATIVE_INT, buffer);
            H5.H5Aclose(attrId);
            return buffer[0];
        }
        catch (HDF5AttributeException e) {
            throw new ObjectNotFoundException("Attribute not found: " + name);
        }
        catch (HDF5Exception ex) {
            String msg = "Error reading attribute from dataset";
            log.error(msg, ex);
            throw new DataAccessException(msg, ex);
        }

    }

    /**
     * Method description
     *
     * @param groupId
     * @param name
     * @return
     */
    public long readLongAttribute(int groupId, String name) {

        try {
            int attrId = H5.H5Aopen_name(groupId, name);
            long[] buffer = new long[1];
            H5.H5Aread(attrId, HDF5Constants.H5T_NATIVE_LONG, buffer);
            H5.H5Aclose(attrId);
            return buffer[0];
        }
        catch (HDF5AttributeException e) {
            throw new ObjectNotFoundException("Attribute not found: " + name);
        }
        catch (HDF5Exception ex) {
            String msg = "Error reading attribute from dataset";
            log.error(msg, ex);
            throw new DataAccessException(msg, ex);
        }

    }

    /**
     * Method description
     *
     * @param groupId
     * @param name
     * @return
     */
    public String readStringAttribute(int groupId, String name) {

        try {
            int attrId = H5.H5Aopen_name(groupId, name);
            int type = H5.H5Aget_type(attrId);
            int size = H5.H5Tget_size(type);
            byte[] buffer = new byte[2 * size];
            H5.H5Aread(attrId, type, buffer);
            H5.H5Aclose(attrId);
            return new String(buffer).trim();
        }
        catch (HDF5AttributeException e) {
            throw new ObjectNotFoundException("Attribute not found: " + name);
        }
        catch (HDF5Exception ex) {
            String msg = "Error reading attribute from dataset";
            log.error(msg, ex);
            throw new DataAccessException(msg, ex);
        }

    }

    /**
     * Method description
     *
     * @param attId
     */
    public void closeAttribute(int attId) {
        try {
            H5.H5Aclose(attId);
        }
        catch (HDF5LibraryException ex) {
            String msg = "Error closing attribute";
            log.error(msg);
            throw new DataAccessException(msg, ex);
        }

    }

    /**
     * Return the dimensions of the dataset.
     */
    private long[] getDims(int datasetId) {
        try {
            int dataspace = H5.H5Dget_space(datasetId);
            int nDims = H5.H5Sget_simple_extent_ndims(dataspace);
            long[] dims = new long[nDims];
            long[] maxDims = new long[nDims];
            H5.H5Sget_simple_extent_dims(dataspace, dims, maxDims);
            return dims;
        }
        catch (HDF5LibraryException ex) {
            String msg = "Error getting dataset dimensions";
            log.error(msg);
            throw new DataAccessException(msg, ex);
        }

    }

    /**
     * Return the dataspace for the dataset
     */
    private int getDataspace(int datasetId) {
        try {
            return H5.H5Dget_space(datasetId);
        }
        catch (HDF5LibraryException ex) {
            String msg = "Error getting dataspace";
            log.error(msg);
            throw new DataAccessException(msg, ex);
        }

    }

    /**
     * Read the entire dataset as an array of strings.
     * for (int i=0; i<rank; i++)
     *
     * @param datasetId
     * @return
     */
    public List<String> readAllStrings(int datasetId) {
        try {
            long[] dims = getDims(datasetId);
            int typeId = H5.H5Dget_type(datasetId);


            int numberOfStrings = 1;
            for (int i = 0; i < dims.length; i++) {
                numberOfStrings *= (int) dims[i];
            }

//          figure out the string size and number of strings
            int stringLength = (int) H5.H5Tget_size(typeId);
            int bufferSize = numberOfStrings * stringLength;
            byte[] byteBuff = new byte[bufferSize];

            // read the string data into byte buff
            int mspace = HDF5Constants.H5S_ALL;
            int fspace = HDF5Constants.H5S_ALL;
            int plist = HDF5Constants.H5P_DEFAULT;
            H5.H5Dread(datasetId, typeId, mspace, fspace, plist, byteBuff);

            // convert byte array into string array
            List<String> strOut = new ArrayList(numberOfStrings);
            for (int i = 0; i < numberOfStrings; i++) {
                strOut.add(new String(byteBuff, i * stringLength, stringLength).trim());
            }

            return strOut;

        }
        catch (HDF5Exception ex) {
            String msg = "Error reading String dataset";
            log.error(msg);
            throw new DataAccessException(msg, ex);
        }

    }

    /**
     * Read the entire dataset as an array of floats.
     *
     * @param datasetId
     * @return
     */
    public float[] readAllFloats(int datasetId) {
        try {
            int size = (int) getDims(datasetId)[0];
            float[] data = new float[size];
            H5.H5Dread(datasetId, HDF5Constants.H5T_NATIVE_FLOAT, HDF5Constants.H5S_ALL,
                    HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, data);
            return data;
        }
        catch (HDF5Exception ex) {
            String msg = "Error reading Float dataset";
            log.error(msg);
            throw new DataAccessException(msg, ex);
        }

    }

    /**
     * Read the entire dataset as an array of ints.  The dsName is
     * relative to the object identfed by groupId (file or group).
     *
     * @param datasetId
     * @return
     */
    public int[] readAllInts(int datasetId) {
        try {
            int size = (int) getDims(datasetId)[0];
            int[] data = new int[size];
            H5.H5Dread(datasetId, HDF5Constants.H5T_NATIVE_INT, HDF5Constants.H5S_ALL,
                    HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, data);
            return data;
        }
        catch (HDF5Exception ex) {
            String msg = "Error reading Int dataset";
            log.error(msg);
            throw new DataAccessException(msg, ex);
        }

    }

    /**
     * Read a section of the dataset as an array of floats.
     *
     * @param datasetId
     * @param fromIndex
     * @param toIndex
     * @return
     */
    public float[] readFloats(int datasetId, int fromIndex, int toIndex) {
        try {

            int dataspace = getDataspace(datasetId);

            int nPts = toIndex - fromIndex + 1;
            long[] offset = new long[1];
            long[] count = new long[1];
            offset[0] = fromIndex;
            count[0] = nPts;
            H5.H5Sselect_hyperslab(dataspace, HDF5Constants.H5S_SELECT_SET, offset, null, count,
                    null);

            // Step 2.   Create a memory dataspace
            long[] memDims = new long[1];
            memDims[0] = nPts;
            int memspace = H5.H5Screate_simple(1, memDims, null);

            // Step 3.  Create a "to" hyperslab
            long[] memOffset = new long[1];
            memOffset[0] = 0;
            long[] memCount = new long[1];
            memCount[0] = nPts;
            H5.H5Sselect_hyperslab(memspace, HDF5Constants.H5S_SELECT_SET, memOffset, null,
                    memCount, null);

            // Step 4.  Read the data
            float[] partialData = new float[nPts];
            H5.H5Dread(datasetId, HDF5Constants.H5T_NATIVE_FLOAT, memspace, dataspace,
                    HDF5Constants.H5P_DEFAULT, partialData);

            return partialData;
        }
        catch (HDF5Exception ex) {
            String msg = "Error reading floats from " + fromIndex + " to " + toIndex;
            log.error(msg);
            throw new DataAccessException(msg, ex);
        }

    }

    /**
     * Read a section of the dataset as an array of floats.
     *
     * @param datasetId
     * @param fromIndex
     * @param toIndex
     * @return
     */
    public double[] readDoubles(int datasetId, int fromIndex, int toIndex) {
        try {

            int dataspace = getDataspace(datasetId);

            int nPts = toIndex - fromIndex + 1;
            long[] offset = new long[1];
            long[] count = new long[1];
            offset[0] = fromIndex;
            count[0] = nPts;
            H5.H5Sselect_hyperslab(dataspace, HDF5Constants.H5S_SELECT_SET, offset, null, count,
                    null);

            // Step 2.   Create a memory dataspace
            long[] memDims = new long[1];
            memDims[0] = nPts;
            int memspace = H5.H5Screate_simple(1, memDims, null);

            // Step 3.  Create a "to" hyperslab
            long[] memOffset = new long[1];
            memOffset[0] = 0;
            long[] memCount = new long[1];
            memCount[0] = nPts;
            H5.H5Sselect_hyperslab(memspace, HDF5Constants.H5S_SELECT_SET, memOffset, null,
                    memCount, null);

            // double  Step 4.  Read the data
            double[] partialData = new double[nPts];
            H5.H5Dread(datasetId, HDF5Constants.H5T_NATIVE_DOUBLE, memspace, dataspace,
                    HDF5Constants.H5P_DEFAULT, partialData);

            return partialData;
        }
        catch (HDF5Exception ex) {
            String msg = "Error reading double from " + fromIndex + " to " + toIndex;
            log.error(msg);
            throw new DataAccessException(msg, ex);
        }

    }

    /**
     * Method description
     *
     * @param datasetId
     * @param fromIndex
     * @param toIndex
     * @return
     */
    public long[] readLongs(int datasetId, int fromIndex, int toIndex) {
        try {

            int dataspace = getDataspace(datasetId);

            int nPts = toIndex - fromIndex + 1;
            long[] offset = new long[1];
            long[] count = new long[1];
            offset[0] = fromIndex;
            count[0] = nPts;
            H5.H5Sselect_hyperslab(dataspace, HDF5Constants.H5S_SELECT_SET, offset, null, count,
                    null);

            // Step 2.   Create a memory dataspace
            long[] memDims = new long[1];
            memDims[0] = nPts;
            int memspace = H5.H5Screate_simple(1, memDims, null);

            // Step 3.  Create a "to" hyperslab
            long[] memOffset = new long[1];
            memOffset[0] = 0;
            long[] memCount = new long[1];
            memCount[0] = nPts;
            H5.H5Sselect_hyperslab(memspace, HDF5Constants.H5S_SELECT_SET, memOffset, null,
                    memCount, null);

            // Step 4.  Read the data
            long[] partialData = new long[nPts];
            H5.H5Dread(datasetId, HDF5Constants.H5T_NATIVE_LLONG, memspace, dataspace,
                    HDF5Constants.H5P_DEFAULT, partialData);

            return partialData;
        }
        catch (HDF5Exception ex) {
            String msg = "Error reading longs from " + fromIndex + " to " + toIndex;
            log.error(msg);
            throw new DataAccessException(msg, ex);
        }

    }

    /**
     * Method description
     *
     * @param dsName
     * @param fromIndex
     * @param toIndex
     * @return
     */
    public int[] readInts(String dsName, int fromIndex, int toIndex) {

        int datasetId = openDataset(dsName);
        int[] ints = readInts(datasetId, fromIndex, toIndex);
        closeDataset(datasetId);
        return ints;


    }

    /**
     * Method description
     *
     * @param datasetId
     * @param fromIndex
     * @param toIndex
     * @return
     */
    public int[] readInts(int datasetId, int fromIndex, int toIndex) {
        try {

            int dataspace = getDataspace(datasetId);

            int nPts = toIndex - fromIndex + 1;
            long[] offset = new long[1];
            long[] count = new long[1];
            offset[0] = fromIndex;
            count[0] = nPts;
            H5.H5Sselect_hyperslab(dataspace, HDF5Constants.H5S_SELECT_SET, offset, null, count,
                    null);

            // Step 2.   Create a memory dataspace
            long[] memDims = new long[1];
            memDims[0] = nPts;
            int memspace = H5.H5Screate_simple(1, memDims, null);

            // Step 3.  Create a "to" hyperslab
            long[] memOffset = new long[1];
            memOffset[0] = 0;
            long[] memCount = new long[1];
            memCount[0] = nPts;
            H5.H5Sselect_hyperslab(memspace, HDF5Constants.H5S_SELECT_SET, memOffset, null,
                    memCount, null);

            // Step 4.  Read the data
            int[] partialData = new int[nPts];
            H5.H5Dread(datasetId, HDF5Constants.H5T_NATIVE_INT, memspace, dataspace,
                    HDF5Constants.H5P_DEFAULT, partialData);

            return partialData;
        }
        catch (HDF5Exception ex) {
            String msg = "Error reading ints from " + fromIndex + " to " + toIndex;
            log.error(msg, ex);
            throw new DataAccessException(msg, ex);
        }

    }


    @Override
    protected void finalize() throws Throwable {
        super.finalize();
        closeFile();
    }
}
