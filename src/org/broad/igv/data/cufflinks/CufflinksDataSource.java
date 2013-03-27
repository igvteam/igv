/*
 * Copyright (c) 2007-2013 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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

package org.broad.igv.data.cufflinks;

import org.apache.commons.math.stat.StatUtils;
import org.broad.igv.Globals;
import org.broad.igv.data.BasicScore;
import org.broad.igv.data.DataSource;
import org.broad.igv.data.GenomeSummaryData;
import org.broad.igv.data.LocusScoreUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.tdf.Accumulator;
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

    public CufflinksDataSource(List<CufflinksValue> valueList, Genome genome) {

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
        return dataMax;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double getDataMin() {
        return dataMin;  //To change body of implemented methods use File | Settings | File Templates.
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
        return TrackType.FPKM;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void setWindowFunction(WindowFunction statType) {
        // Ignored
    }

    @Override
    public boolean isLogNormalized() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public WindowFunction getWindowFunction() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
