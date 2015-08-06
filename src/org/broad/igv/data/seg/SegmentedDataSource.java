/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
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
