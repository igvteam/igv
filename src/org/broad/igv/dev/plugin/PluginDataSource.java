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

package org.broad.igv.dev.plugin;

import org.apache.log4j.Logger;
import org.broad.igv.data.DataSource;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.DataTrack;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;

import java.io.IOException;
import java.util.*;

/**
 * A feature source which derives its information
 * from a command line plugin
 * User: jacob
 * Date: 2012/05/01
 */
public class PluginDataSource extends PluginSource<LocusScore, LocusScore> implements DataSource {

    private static Logger log = Logger.getLogger(PluginDataSource.class);

    private double dataMin;
    private double dataMax;
    private WindowFunction windowFunction;

    public PluginDataSource(List<String> commands, LinkedHashMap<Argument, Object> arguments, Map<String, String> parsingAttrs, String specPath) {
        super(commands, arguments, parsingAttrs, specPath);
    }

    protected String createTempFile(Track track, Argument argument, String chr, int start, int end, int zoom) throws IOException {
        DataTrack dataTrack = (DataTrack) track;
        List<LocusScore> features = dataTrack.getSummaryScores(chr, start, end, zoom);
        return super.createTempFile(features, argument);
    }

    @Override
    public List<LocusScore> getSummaryScoresForRange(String chr, int startLocation, int endLocation, int zoom) {
        try {
            Iterator<LocusScore> iter = super.getFeatures(chr, startLocation, endLocation, zoom);
            List<LocusScore> list = new ArrayList<LocusScore>(1000);
            LocusScore score;
            dataMin = Double.MAX_VALUE;
            dataMax = -Double.MAX_VALUE;
            while (iter.hasNext()) {
                score = iter.next();
                dataMin = Math.min(dataMin, score.getScore());
                dataMax = Math.max(dataMax, score.getScore());
                list.add(score);
            }
            return list;
        } catch (IOException e) {
            log.error(e);
            throw new RuntimeException(e);
        }
    }


    @Override
    public double getDataMax() {
        return dataMax;
    }

    @Override
    public double getDataMin() {
        return dataMin;
    }

    @Override
    public TrackType getTrackType() {
        return TrackType.PLUGIN;
    }

    @Override
    public void setWindowFunction(WindowFunction statType) {
        this.windowFunction = statType;
    }

    @Override
    public WindowFunction getWindowFunction() {
        return windowFunction;
    }

    @Override
    public boolean isLogNormalized() {
        return false; //TODO
    }

    @Override
    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return null; //TODO
    }
}
