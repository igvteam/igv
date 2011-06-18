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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.h5;

import ncsa.hdf.hdf5lib.H5;
import ncsa.hdf.hdf5lib.exceptions.HDF5LibraryException;
import ncsa.hdf.hdf5lib.exceptions.HDF5SymbolTableException;
import org.apache.log4j.Logger;

import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 */
public class HDF5EntityCache {

    private Logger log = Logger.getLogger(HDF5EntityCache.class);

    int fileId;

    Map<String, Integer> datasetCache;

    Map<String, Integer> groupCache;

    public HDF5EntityCache(int fileId) {
        this.fileId = fileId;
        datasetCache = new HashMap();
        groupCache = new HashMap();
    }

    public synchronized void closeAllEntities() {
        for (Integer dsId : datasetCache.values()) {
            closeDataset(dsId);
        }

        datasetCache.clear();

    }

    public int openDataset(String dsName) {
        try {
            Integer dsId = datasetCache.get(dsName);
            if (dsId == null) {
                dsId = H5.H5Dopen(fileId, dsName);
                datasetCache.put(dsName, dsId);
            }
            return dsId;
        } catch (HDF5SymbolTableException ex) {
            throw new ObjectNotFoundException("Error opening dataset: " + dsName);
        } catch (HDF5LibraryException ex) {
            log.error("Error opening dataset", ex);
            throw new RuntimeException(ex);
        }

    }

    /**
     * Close an hdf5 dataset
     */
    private void closeDataset(int datasetId) {

        try {
            H5.H5Dclose(datasetId);
        } catch (HDF5LibraryException ex) {
            log.error("Error closing dataset", ex);
            throw new RuntimeException("Error closing dataset");
        }

    }


    @Override
    protected void finalize() throws Throwable {
        try {
            closeAllEntities();
        } finally {
            super.finalize();
        }
    }
    /*
    static class DatasetReference {
    int datasetId;
    DatasetReference(int datasetId) {
    datasetId = datasetId;
    }
    }
     * */
}
