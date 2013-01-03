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


package org.broad.igv.track;


import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.data.CoverageDataSource;
import org.broad.igv.data.DataSource;
import org.broad.igv.feature.LocusScore;
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

    public DataSourceTrack(ResourceLocator locator, String id, String name, DataSource dataSource) {
        super(locator, id, name);

        this.dataSource = dataSource;
        setTrackType(dataSource.getTrackType());
        List<LocusScore> scores = this.dataSource.getSummaryScoresForRange(Globals.CHR_ALL, -1, -1, -1);

        float min = (float) dataSource.getDataMin();
        float max = (float) dataSource.getDataMax();
        float baseline = 0;

        // If the range is all + numbers set the min to zero
        if (min > 0) {
            min = 0;
        }
        for (LocusScore score : scores) {
            max = Math.max(max, score.getScore());
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
