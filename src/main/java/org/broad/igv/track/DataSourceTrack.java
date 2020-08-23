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


package org.broad.igv.track;


import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.data.CombinedDataSource;
import org.broad.igv.data.CoverageDataSource;
import org.broad.igv.data.DataSource;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.util.ResourceLocator;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.Collection;
import java.util.Collections;
import java.util.List;


/**
 * @author jrobinso
 */

public class DataSourceTrack extends DataTrack {

    private static Logger log = Logger.getLogger(DataSourceTrack.class);

    public DataSource dataSource;


    public DataSourceTrack(ResourceLocator locator, String id, String name, DataSource dataSource) {
        super(locator, id, name);

        this.dataSource = dataSource;

        if (this.dataSource != null) {
            setTrackType(dataSource.getTrackType());
            List<LocusScore> scores = this.dataSource.getSummaryScoresForRange(Globals.CHR_ALL, -1, -1, 0);

            if (scores.size() > 0) {
                initScale(dataSource, scores);
            }
        }
    }

    public DataSourceTrack() {

    }


    void initScale(DataSource dataSource, List<LocusScore> scores) {

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

    public LoadedDataInterval<List<LocusScore>> getSummaryScores(String chr, int startLocation, int endLocation, int zoom) {

        List<LocusScore> tmp = dataSource.getSummaryScoresForRange(chr, startLocation, endLocation, zoom);
        if (tmp == null) tmp = Collections.EMPTY_LIST;
        if (dataRange == null) {
            initScale(dataSource, tmp);
        }
        return new LoadedDataInterval<>(chr, startLocation, endLocation, zoom, tmp);
    }


    @Override
    public void setWindowFunction(WindowFunction statType) {
        clearCaches();
        if (dataSource != null) {
            dataSource.setWindowFunction(statType);
        }
    }

    public boolean isLogNormalized() {
        return dataSource.isLogNormalized();
    }


    public WindowFunction getWindowFunction() {

        return dataSource != null ? dataSource.getWindowFunction() : null;
    }


    @Override
    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return dataSource != null ? dataSource.getAvailableWindowFunctions() : null;
    }

    public void updateTrackReferences(List<Track> allTracks) {
        if (dataSource instanceof CombinedDataSource) {
            ((CombinedDataSource) dataSource).updateTrackReferences(allTracks);
        }
    }

    @Override
    public void dispose() {
        super.dispose();
        if (dataSource != null) {
            dataSource.dispose();
            dataSource = null;
        }
    }


    public void marshalXML(Document document, Element element) {

        super.marshalXML(document, element);

        if (dataSource != null) {
            if (dataSource instanceof CoverageDataSource) {
                boolean normalize = ((CoverageDataSource) dataSource).getNormalize();
                if (normalize) {
                    element.setAttribute("normalize", "true");
                }
            }
        }
    }

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);

        if (dataSource != null) {
            if(dataSource instanceof CoverageDataSource && element.hasAttribute("normalize")){
                ((CoverageDataSource) dataSource).setNormalize(Boolean.parseBoolean(element.getAttribute("normalize")));
            }
        }


    }
}
