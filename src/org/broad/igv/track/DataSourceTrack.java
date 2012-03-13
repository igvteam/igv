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


package org.broad.igv.track;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.data.CoverageDataSource;
import org.broad.igv.data.DataSource;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.util.ResourceLocator;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;


/**
 * @author jrobinso
 */
public class DataSourceTrack extends DataTrack {

    private static Logger log = Logger.getLogger(DataSourceTrack.class);

    private DataSource dataSource;
    boolean normalize = false;

    public DataSourceTrack(ResourceLocator locator, String id, String name, DataSource dataSource, Genome genome) {
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

    public List<LocusScore> getSummaryScores(String chr, int startLocation, int endLocation, int zoom) {
        List<LocusScore> tmp = dataSource.getSummaryScoresForRange(chr, startLocation, endLocation, zoom);
        return tmp == null ? new ArrayList() : tmp;

    }


    @Override
    public void setWindowFunction(WindowFunction statType) {
        clearCaches();
        dataSource.setWindowFunction(statType);

    }

    public boolean isLogNormalized() {
        return dataSource.isLogNormalized();
    }


    public WindowFunction getWindowFunction() {
        return dataSource.getWindowFunction();
    }


    @Override
    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return dataSource.getAvailableWindowFunctions();
    }

    @Override
    public Map<String, String> getPersistentState() {
        Map<String, String> properties = super.getPersistentState();
        if (normalize != false) {
            properties.put("normalize", String.valueOf(normalize));
        }
        return properties;
    }


    @Override
    public void restorePersistentState(Map<String, String> attributes) {
        super.restorePersistentState(attributes);
        String as = attributes.get("normalize");
        if (as != null) {
            try {
                normalize = Boolean.parseBoolean(as);
                if (dataSource != null && dataSource instanceof CoverageDataSource) {
                    ((CoverageDataSource) dataSource).setNormalize(normalize);
                }
            } catch (Exception e) {
                log.error("Error restoring session.  Invalid normalization value: " + normalize);
            }
        }
    }
}
