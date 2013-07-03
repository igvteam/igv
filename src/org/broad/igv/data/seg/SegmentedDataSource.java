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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.data.seg;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.Globals;
import org.broad.igv.data.DataSource;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 *
 * @author jrobinso
 */
public class SegmentedDataSource implements DataSource {

    /**
     * Identifies the track.  Often a sample name.
     */
    String trackIdentifier;
    SegmentedDataSet dataset;

    /**
     * Constructs ...
     *
     * @param trackIdentifier
     * @param dataset
     */
    public SegmentedDataSource(String trackIdentifier, SegmentedDataSet dataset) {
        this.trackIdentifier = trackIdentifier;
        this.dataset = dataset;
    }

    public TrackType getTrackType() {
        return dataset.getType();
    }


    public void setWindowFunction(WindowFunction statType) {

        // ignore
    }


    public boolean isLogNormalized() {
        return dataset.isLogNormalized();
    }


    public double getDataMax() {
        return dataset.getDataMax(Globals.CHR_ALL);
    }


    public double getDataMin() {
        return dataset.getDataMin(Globals.CHR_ALL);
    }


    public double getMedian(int zoom, String chr) {
        return 1.0;
    }

    private List<LocusScore> getSegments(String chr) {

        return dataset.getSegments(trackIdentifier, chr);

    }


    public List<LocusScore> getSummaryScoresForRange(String chr, int startLocation,
                                                     int endLocation, int zoom) {
        if (chr.equals(Globals.CHR_ALL)) {
            return getWholeGenomeScores();
        }
        return getSegments(chr);
    }


    public List<LocusScore> getWholeGenomeScores() {
        return dataset.getWholeGenomeScores(trackIdentifier);
    }


    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return new ArrayList();
    }

    @Override
    public void dispose() {

    }


    public WindowFunction getWindowFunction() {
        return null;
    }

    public String getTrackIdentifier() {
        return trackIdentifier;
    }
}
