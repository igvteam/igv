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

package org.broad.igv.h5;

import ncsa.hdf.hdf5lib.H5;
import ncsa.hdf.hdf5lib.HDF5Constants;

/**
 * Test of reading HDF5 string.
 */
public class TestH5String {
    public static void main(String[] argv) {
        int fid = -1, did = -1, tid = -1, sid = -1;
        long[] dims = {2};
        int rank = 1;
        int strLength = 14;
        String[] strIn = {"test string 1 ", "test string 2 "};
        String[] strOut = null;
        byte[] byteBuff = null;

        try {
            // create a new file
            fid = H5.H5Fcreate("testH5string.h5",
                    HDF5Constants.H5F_ACC_TRUNC,
                    HDF5Constants.H5P_DEFAULT,
                    HDF5Constants.H5P_DEFAULT);

            // create a new dataset of two strings
            tid = H5.H5Tcopy(HDF5Constants.H5T_C_S1);
            H5.H5Tset_size(tid, strLength);
            sid = H5.H5Screate_simple(1, dims, null);
            did = H5.H5Dcreate(fid, "/string", tid, sid, HDF5Constants.H5P_DEFAULT);
            byteBuff = (strIn[0] + strIn[1]).getBytes(); // for large number of strings, use StringBuffer
            H5.H5Dwrite(did, tid, HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, byteBuff);

            try {
                H5.H5Tclose(tid);
            } catch (Exception ex) {
            }
            try {
                H5.H5Sclose(sid);
            } catch (Exception ex) {
            }
            try {
                H5.H5Dclose(did);
            } catch (Exception ex) {
            }


            // read the string from file and print it out.
            did = H5.H5Dopen(fid, "/string");
            tid = H5.H5Dget_type(did);
            sid = H5.H5Dget_space(did);

            int numberOfStrings = 1;
            for (int i = 0; i < rank; i++)
                numberOfStrings *= (int) dims[i];

            // figure out the string size and number of strings
            int stringLength = (int) H5.H5Tget_size(tid);
            int bufferSize = numberOfStrings * stringLength;
            byteBuff = new byte[bufferSize];

            // read the string data into byte buff
            int mspace = HDF5Constants.H5S_ALL;
            int fspace = HDF5Constants.H5S_ALL;
            int plist = HDF5Constants.H5P_DEFAULT;
            H5.H5Dread(did, tid, mspace, fspace, plist, byteBuff);

            // convert byte array into string array
            strOut = new String[numberOfStrings];
            for (int i = 0; i < numberOfStrings; i++) {
                strOut[i] = new String(byteBuff, i * stringLength, stringLength).trim();
            }
        } catch (Exception ex) {
            System.out.println(ex);
        } finally {
            try {
                H5.H5Tclose(tid);
            } catch (Exception ex) {
            }
            try {
                H5.H5Sclose(sid);
            } catch (Exception ex) {
            }
            try {
                H5.H5Dclose(did);
            } catch (Exception ex) {
            }
            try {
                H5.H5Fclose(fid);
            } catch (Exception ex) {
            }
        }

        //print out string array
        if (strOut != null) {
            for (int i = 0; i < strOut.length; i++)
                System.out.println(strOut[i]);
        }
    }
}

