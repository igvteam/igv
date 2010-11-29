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
 * FeatureTrackH5.java
 *
 * Created on November 12, 2007, 8:22 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package org.broad.igv.track;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.data.DataSource;
import org.broad.igv.data.HDFDataManager;
import org.broad.igv.data.HDFDataSource;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.h5.ObjectNotFoundException;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Represents a feature track backed by an HDF-5 file.
 *
 * @author jrobinso
 */
public class HDFDataTrack extends DataTrack {

    private static Logger log = Logger.getLogger(HDFDataTrack.class);
    private DataSource dataSource;
    int trackNumber;
    private static Color K4_COLOR = new Color(0, 150, 0);
    private static Color K9_COLOR = new Color(100, 0, 0);
    private static Color K27_COLOR = Color.RED;
    private static Color K36_COLOR = new Color(0, 0, 150);

    /**
     * @param dataManager
     * @param locator
     * @param name
     * @param trackNumber
     * @throws FileNotFoundException
     */
    public HDFDataTrack(HDFDataManager dataManager,
                        ResourceLocator locator,
                        String name,
                        int trackNumber)
            throws FileNotFoundException {

        super(locator, name + trackNumber, name);
        //this.windowFunction = scoreType;
        this.trackNumber = trackNumber;

        dataSource = new HDFDataSource(dataManager, name, trackNumber);


        // By default autoscale so that (1) lower limit is zero,  or (2) scale
        // is symetrical

        float max = (float) dataSource.getDataMax();
        float min = Math.min(0, (float) dataSource.getDataMin());
        if (min < 0) {
            float absMax = Math.max(Math.abs(max), Math.abs(min));
            min = -absMax;
            max = absMax;
        }
        if (max == min) {
            max += 1;
        }

        float baseline = 0;

        setDataRange(new DataRange(min, baseline, max));

        try {
            setTrackType(dataSource.getTrackType());

        } catch (Exception exception) {

            // Just log the error.  This is not a fatal exception, the
            // track will be set to the "GNERIC" type
            log.error("Unknown track type: " + dataSource.getTrackType());
        }

        // Temporary color hack until attribute file supports color

        String tmp = getName();
        if (tmp.contains("K4")) {
            setColor(K4_COLOR);
        } else if (tmp.contains("K9")) {
            setColor(K9_COLOR);
        } else if (tmp.contains("K27")) {
            setColor(K27_COLOR);
        } else if (tmp.contains("K36")) {
            setColor(K36_COLOR);
        }
    }

    // TODO add zoom as an argument, remove reference to IGVMainFrame
    /**
     * Method description
     *
     * @param chr
     * @param startLocation
     * @param endLocation
     * @param zoom
     * @return
     */
    public List<LocusScore> getSummaryScores(String chr, int startLocation,
                                             int endLocation, int zoom) {
        try {
            return dataSource.getSummaryScoresForRange(chr, startLocation, endLocation, zoom);

        } catch (ObjectNotFoundException ex) {

            // This is an expected condition,  return empty list
            return new ArrayList();
        }

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
     * @param statType
     */
    public void setStatType(WindowFunction statType) {
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

    //  TODO -- generalize.  Move calculation into pluggable class
    /**
     * Method description
     *
     * @param chr
     * @param start
     * @param end
     * @param zoom
     * @param type
     * @return
     */
    public float getRegionScore(String chr, int start, int end, int zoom,
                                RegionScoreType type) {

        if (end <= start) {
            return 0;
        }

        if (isRegionScoreType(type)) {

            List<LocusScore> scores = getSummaryScores(chr, start, end, zoom);
            float regionScore = 0;
            int intervalSum = 0;
            for (LocusScore score : scores) {
                if ((score.getEnd() >= start) && (score.getStart() <= end)) {
                    int interval = Math.min(end, score.getEnd()) - Math.max(start, score.getStart());
                    float value = score.getScore();
                    regionScore += value * interval;
                    intervalSum += interval;
                }

            }
            if (intervalSum <= 0) {
                return -Float.MAX_VALUE;
            } else {
                regionScore /= intervalSum;
                return (type == RegionScoreType.DELETION)
                        ? -regionScore : regionScore;
            }

        } else {
            return -Float.MAX_VALUE;
        }

    }

    /**
     * Refresh the underlying data for the track.
     *
     * @param timestamp used to prevent multiple refreshes from the same request
     */
    public void refreshData(long timestamp) {
        dataSource.refreshData(timestamp);

    }

    /**
     * Method description
     *
     * @param chr
     * @param startLocation
     * @param endLocation
     * @return
     */
    public List<List<IGVFeature>> getFeaturesByLevels(String chr,
                                                   int startLocation, int endLocation) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * Method description
     *
     * @return
     */
    public int getNumberOfFeatureLevels() {
        return 1;
    }

    /**
     * Method description
     *
     * @param isMultiLevel
     */
    public void setMultiLevelFeatures(boolean isMultiLevel) {
        // ignore
    }

    /**
     * Method description
     *
     * @return
     */
    public boolean isMuliLevelFeatures() {
        return false;
    }

    // TODO -- get window functions dynamically from data
    static List<WindowFunction> wfs = new ArrayList();

    static {
        wfs.add(WindowFunction.percentile10);
        wfs.add(WindowFunction.median);
        wfs.add(WindowFunction.mean);
        wfs.add(WindowFunction.percentile90);
        wfs.add(WindowFunction.max);

    }

    @Override
    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return wfs;
    }
}
