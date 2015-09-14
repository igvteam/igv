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

package org.broad.igv.data.cufflinks;

import org.apache.commons.math.stat.StatUtils;
import org.broad.igv.Globals;
import org.broad.igv.data.BasicScore;
import org.broad.igv.data.DataSource;
import org.broad.igv.data.GenomeSummaryData;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.util.collections.DownsampledDoubleArrayList;

import java.util.*;

/**
 * @author jrobinso
 *         Date: 3/8/13
 *         Time: 2:52 PM
 */
public class CufflinksDataSource implements DataSource {

    double dataMax;
    double dataMin;
    Map<String, List<LocusScore>> values;
    Map<String, String> chrAliasMap;
    List<LocusScore> wholeGenomeScores;

    private int sampleIndex = 0;

    /** Constructor for files with multiple samples. Cufflinks files have multiple
     * columns at end (just FPKM files as of this writing) representing each sample
     * @param sampleIndex
     * @param allValues
     * @param genome
     */
    public CufflinksDataSource(int sampleIndex, List<FPKMValue> allValues, Genome genome){
        this(getSampleValues(sampleIndex, allValues), genome);
        this.sampleIndex = sampleIndex;
    }

    private static List<? extends LocusScore> getSampleValues(int sampleIndex, List<FPKMValue> allValues) {
        List<FPKMSampleValue> sampleValueList = new ArrayList<FPKMSampleValue>(allValues.size());
        for(FPKMValue value: allValues){
            sampleValueList.add(value.getSampleValue(sampleIndex));
        }
        return sampleValueList;
    }

    public CufflinksDataSource(List<? extends LocusScore> valueList, Genome genome) {

        chrAliasMap = new HashMap<String, String>();
        values = new HashMap<String, List<LocusScore>>();

        DownsampledDoubleArrayList sampledData = sampleValues(valueList, genome);

        // Sort
        for (List<LocusScore> chrValues : values.values()) {
            FeatureUtils.sortFeatureList(chrValues);
        }

        double[] sd = sampledData.toArray();
        if (sd.length > 0) {
            dataMin = Math.min(0, StatUtils.percentile(sd, 5));
            dataMax = StatUtils.percentile(sd, 95);
        } else {
            dataMin = 0;
            dataMax = 100;
        }

        calculateWholeGenomeScores(genome);
    }

    private void calculateWholeGenomeScores(Genome genome){
        GenomeSummaryData genomeSummaryData = new GenomeSummaryData(genome, new String[]{"*"});
        for (Map.Entry<String, List<LocusScore>> entry : values.entrySet()) {
            String chr = entry.getKey();
            List<LocusScore> scores = entry.getValue();
            int[] positions = new int[scores.size()];
            float[] values = new float[scores.size()];
            for (int i = 0; i < scores.size(); i++) {
                LocusScore s = scores.get(i);
                positions[i] = s.getStart();
                values[i] = s.getScore();
            }
            Map<String, float[]> tmp = new HashMap<String, float[]>(1);
            tmp.put("*", values);
            genomeSummaryData.addData(chr, positions, tmp);
        }
        int[] positions = genomeSummaryData.getLocations();
        float[] values = genomeSummaryData.getData("*");
        wholeGenomeScores = new ArrayList<LocusScore>(positions.length);
        for (int i = 0; i < positions.length; i++) {
            wholeGenomeScores.add(new BasicScore(positions[i], positions[i] + 1, values[i]));
        }

    }

    /**
     * Sample the first 10,000 values to set scale
     * Also separate data into chromosomes
     * @param valueList
     * @param genome
     */
    private DownsampledDoubleArrayList sampleValues(List<? extends LocusScore> valueList, Genome genome){
        DownsampledDoubleArrayList sampledData = new DownsampledDoubleArrayList(5000, 10000);
        for (LocusScore val : valueList) {
            String chr = val.getChr();

            List<LocusScore> chrValues = values.get(chr);
            if (chrValues == null) {
                chrValues = new ArrayList<LocusScore>();
                values.put(chr, chrValues);
                if (genome != null) {
                    String alias = genome.getChromosomeAlias(chr);
                    chrAliasMap.put(alias, chr);
                }


            }
            sampledData.add(val.getScore());
            chrValues.add(val);
        }
        return sampledData;
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
    public List<LocusScore> getSummaryScoresForRange(String chr, int startLocation, int endLocation, int zoom) {

        if (chr.equals(Globals.CHR_ALL)) {
            return wholeGenomeScores;
        }

        if (chrAliasMap.containsKey(chr)) {
            return values.get(chrAliasMap.get(chr));
        } else {
            return values.get(chr);
        }
    }

    @Override
    public TrackType getTrackType() {
        return TrackType.FPKM;
    }

    @Override
    public void setWindowFunction(WindowFunction statType) {
        // Ignored
    }

    @Override
    public boolean isLogNormalized() {
        return false;
    }

    @Override
    public WindowFunction getWindowFunction() {
        return null;
    }

    @Override
    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return Arrays.asList(WindowFunction.none);
    }

    @Override
    public void dispose() {

    }
}
