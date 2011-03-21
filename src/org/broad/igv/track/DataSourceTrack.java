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


package org.broad.igv.track;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.data.DataSource;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.h5.ObjectNotFoundException;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.ui.IGVMainFrame;
import org.broad.igv.util.ResourceLocator;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;


/**
 * @author jrobinso
 */
public class DataSourceTrack extends DataTrack {

    private DataSource dataSource;
    public static double log2 = Math.log(2);

    // private WindowFunction windowFunction = WindowFunction.median;

    /**
     * Constructs ...
     *
     * @param locator
     * @param name
     * @param dataSource
     */
    public DataSourceTrack(ResourceLocator locator, String id, String name, DataSource dataSource) {
        super(locator, id, name);

        this.dataSource = dataSource;
        setTrackType(dataSource.getTrackType());

        float min = (float) dataSource.getDataMin();
        float max = (float) dataSource.getDataMax();
        float baseline = 0;
        if (min > 0) {
            min = 0;
        }

        setDataRange(new DataRange(min, baseline, max));

    }

    public DataSource getDataSource() {
        return dataSource;
    }

    /**
     * Method description
     *
     * @param chr
     * @param startLocation
     * @param endLocation
     * @param zoom
     * @return
     */
    public List<LocusScore> getSummaryScores(String chr, int startLocation, int endLocation, int zoom) {
        try {
            List<LocusScore> tmp = dataSource.getSummaryScoresForRange(chr, startLocation, endLocation, zoom);
            return tmp == null ? new ArrayList() : tmp;
        }
        catch (ObjectNotFoundException ex) {
            return new ArrayList();
        }
    }

    /**
     * Method description
     *
     * @param statType
     */
    @Override
    public void setWindowFunction(WindowFunction statType) {
        clearCaches();
        dataSource.setWindowFunction(statType);

    }


    /**
     * Method description
     *
     * @return
     */
    public boolean isLogNormalized() {
        return dataSource.isLogNormalized();
    }

    /**
     * Method description
     *
     * @param timestamp
     */
    @Override
    public void refreshData(long timestamp) {
        dataSource.refreshData(timestamp);
    }

    /**
     * Method description
     *
     * @return
     */
    public WindowFunction getWindowFunction() {
        return dataSource.getWindowFunction();
    }

    /**
     * Method description
     *
     * @return
     */
    @Override
    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return dataSource.getAvailableWindowFunctions();
    }


}
