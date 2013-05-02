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

package org.broad.igv.data.cufflinks;

import org.apache.commons.math.stat.StatUtils;
import org.broad.igv.Globals;
import org.broad.igv.data.BasicScore;
import org.broad.igv.data.DataSource;
import org.broad.igv.data.GenomeSummaryData;
import org.broad.igv.data.LocusScoreUtils;
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

    public CufflinksDataSource(List<? extends CufflinksValue> valueList, Genome genome) {

        chrAliasMap = new HashMap<String, String>();
        // Sample the first 10,000 values to set scale
        DownsampledDoubleArrayList sampledData = new DownsampledDoubleArrayList(5000, 10000);
        values = new HashMap<String, List<LocusScore>>();

        for (CufflinksValue val : valueList) {
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

        // Sort
        for (List<LocusScore> chrValues : values.values()) {
            LocusScoreUtils.sortFeatureList(chrValues);
        }

        double[] sd = sampledData.toArray();
        if (sd.length > 0) {
            dataMin = Math.min(0, StatUtils.percentile(sd, 5));
            dataMax = StatUtils.percentile(sd, 95);
        } else {
            dataMin = 0;
            dataMax = 100;
        }

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
}
