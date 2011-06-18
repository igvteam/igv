/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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


package org.broad.igv.synteny;


import org.broad.igv.data.DataUtils;
import org.broad.igv.data.WiggleDataset;
import org.broad.igv.data.WiggleParser;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.util.ResourceLocator;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;


/**
 * @author jrobinso
 */
public class BlastMappingUtils {


    /**
     * Method description
     *
     * @param dataset
     * @param mappings
     * @param outputFile
     */
    public static void MapWigFile(WiggleDataset dataset, List<BlastMapping> mappings, File outputFile) {

        PrintWriter pw = null;

        try {
            pw = new PrintWriter(new FileWriter(outputFile));

            // Segregate and sort ortholog
            Map<String, List<BlastMapping>> featureMap = new LinkedHashMap();
            for (BlastMapping f : mappings) {
                List<BlastMapping> fList = featureMap.get(f.getChr());
                if (fList == null) {
                    fList = new ArrayList();
                    featureMap.put(f.getChr(), fList);
                }
                fList.add(f);
            }

            for (List<BlastMapping> featureList : featureMap.values()) {
                FeatureUtils.sortFeatureList(featureList);
            }

            pw.println("Chr\tStart\tEnd\t\tnucCount (Kwal)");
            // Loop through chromosomes
            for (String chr : featureMap.keySet()) {
                List<BlastMapping> mappingList = featureMap.get(chr);
                List<BasicScore> scores = new ArrayList(mappingList.size());

                for (BlastMapping mapping : mappingList) {
                    BlastMapping.Block queryBlock = mapping.getQueryBlock();
                    BlastMapping.Block subjectBlock = mapping.getSubjectBlock();

                    int[] subjectPosition = dataset.getStartLocations(subjectBlock.getContig());
                    float[] data = dataset.getData("ignore", subjectBlock.getContig());

                    int s0 = subjectBlock.getStart();
                    int s1 = subjectBlock.getEnd();
                    int q0 = queryBlock.getStart();
                    int q1 = queryBlock.getEnd();

                    int idx0 = DataUtils.getIndexBefore(subjectPosition, Math.min(s0, s1));
                    int idx1 = DataUtils.getIndexBefore(subjectPosition, Math.max(s0, s1)) + 1;
                    if (idx1 == subjectPosition.length) {
                        idx1--;
                    }

                    double beta = ((double) (q1 - q0)) / (s1 - s0);

                    for (int i = idx0; i <= idx1; i++) {
                        int pos = (int) (q0 + beta * (subjectPosition[i] - s0));
                        float d = data[i];
                        scores.add(new BasicScore(chr, pos, pos + 1, d));
                    }
                }

                FeatureUtils.sortFeatureList(scores);
                for (BasicScore s : scores) {
                    pw.println(s.getChr() + "\t" + s.getStart() + "\t" +
                            (s.getEnd() + 1) + "\t\t" + s.getScore());

                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        finally {
            pw.close();
        }
    }


    static class BasicScore implements LocusScore {

        String chr;
    int start;
    int end;
    float score;

    public BasicScore(String chromosome, int start, int end, float score) {
        this.chr = chromosome;
        this.start = start;
        this.end = end;
        this.score = score;
    }

    public BasicScore(BasicScore bs) {
        this.chr = bs.chr;
        this.start = bs.start;
        this.end = bs.end;
        this.score = bs.score;
    }

    public BasicScore copy() {
        return new BasicScore(this);
    }

    public String getChr() {
        return chr;
    }



    public int getStart() {
        return start;
    }

    public float getScore() {
        return score;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public String getValueString(double position, WindowFunction windowFunction) {
        return "Value: " + score;
        //throw new UnsupportedOperationException("Not supported yet.");
    }

    public int getExtendedStart() {
        return getStart();
    }

    public int getExtendedEnd() {
        return getEnd();
    }
    }


}
