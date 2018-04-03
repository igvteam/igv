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

package org.broad.igv.cli_plugin;

import org.apache.log4j.Logger;
import org.broad.igv.data.AbstractDataSource;
import org.broad.igv.data.DataTile;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;
import htsjdk.tribble.Feature;

import java.io.IOException;
import java.util.*;

/**
 * A data source which derives its information
 * from a command line cli_plugin. Note that the input
 * User: jacob
 * Date: 2012/05/01
 */
public class PluginDataSource extends AbstractDataSource {

    private static Logger log = Logger.getLogger(PluginDataSource.class);
    private double dataMin;
    private double dataMax;

    private PluginSource<Feature, LocusScore> pluginFeatureSource;

    private Map<String, Integer> longestFeatureMap = new HashMap<String, Integer>();

    public PluginDataSource(Genome genome, List<String> commands, LinkedHashMap<Argument, Object> arguments, PluginSpecReader.Output outputAttr, String specPath) {
        super(genome);
        pluginFeatureSource = new PluginFeatureSource(commands, arguments, outputAttr, specPath);
    }

//    protected String createTempFile(Track track, Argument argument, String chr, int start, int end, int zoom) throws IOException {
//        if(track instanceof DataTrack){
//            DataTrack dataTrack = (DataTrack) track;
//            List<LocusScore> features = dataTrack.getSummaryScores(chr, start, end, zoom);
//            return pluginFeatureSource.createTempFile(features, argument);
//        } else if (track instanceof AlignmentTrack) {
//            List<Alignment> alignments = pluginFeatureSource.getAlignmentsForRange((AlignmentTrack) track, chr, start, end, zoom);
//            return pluginFeatureSource.createTempFile(alignments, argument);
//        }
//        throw new IllegalArgumentException("Unsupported track type: " + track.getClass());
//    }

    @Override
    protected DataTile getRawData(String chr, int startLocation, int endLocation) {

        int queryLength = endLocation - startLocation;
        int longestFeature = getLongestFeature(chr);
        int adjustedStart = startLocation;
        int adjustedEnd = endLocation;
        if(queryLength < longestFeature){
            int halfDiff = (longestFeature - queryLength)/2;
            adjustedStart = startLocation - halfDiff - 1;
            adjustedEnd = endLocation + halfDiff;
        }

        try {
            Iterator<LocusScore> iter = pluginFeatureSource.getFeatures(chr, adjustedStart, adjustedEnd, -1);
            List<LocusScore> list = new ArrayList<LocusScore>(1000);
            LocusScore score;
            dataMin = Double.MAX_VALUE;
            dataMax = -Double.MAX_VALUE;

            while (iter.hasNext()) {
                score = iter.next();
                dataMin = Math.min(dataMin, score.getScore());
                dataMax = Math.max(dataMax, score.getScore());
                list.add(score);
                longestFeature = Math.max(longestFeature, score.getEnd() - score.getStart());
            }
            longestFeatureMap.put(chr, longestFeature);

            int length = list.size();
            int[] startLocations = new int[length];
            int[] endLocations = new int[length] ;
            float[] scores = new float[length];
            String[] names = new String[length];
            int idx = 0;
            for(LocusScore locusScore: list){
                startLocations[idx] = locusScore.getStart();
                endLocations[idx] = locusScore.getEnd();
                scores[idx] = locusScore.getScore();
                names[idx] = Locus.getFormattedLocusString(locusScore.getChr(), locusScore.getStart(), locusScore.getEnd());
                idx++;
            }
            return new DataTile(startLocations, endLocations, scores, names);

        } catch (IOException e) {
            log.error(e.getMessage(), e);
            throw new RuntimeException(e);
        }
    }

    @Override
    protected List<LocusScore> getPrecomputedSummaryScores(String chr, int startLocation, int endLocation, int zoom) {
        return null;
    }

    @Override
    public int getLongestFeature(String chr) {
        return longestFeatureMap.containsKey(chr) ? longestFeatureMap.get(chr) : 1000;
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
    public boolean isLogNormalized() {
        return false; //TODO
    }

    @Override
    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return Arrays.asList(WindowFunction.none); //TODO
    }

//    public void setQueryTracker(QueryTracker queryTracker) {
//        this.pluginFeatureSource.setQueryTracker(queryTracker);
//    }
}
