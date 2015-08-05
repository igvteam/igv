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

package org.broad.igv.peaks;

import org.broad.igv.Globals;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.WindowFunction;
import htsjdk.tribble.Feature;

/**
 * Note:  implementing tribble.Feature will allow us to index these files in the future.
 *
 * @author jrobinso
 * @date Apr 22, 2011
 */
public class Peak implements LocusScore, htsjdk.tribble.Feature {

    String chr;
    int start;
    int end;
    private String name;
    private float combinedScore;
    private float[] timeScores;
    boolean dynamic = false;
    float dynamicScore;
    private float foldChange;

    public Peak(String chr, int start, int end, String name, float combinedScore, float[] timeScores) {
        this.chr = chr;
        this.combinedScore = combinedScore;
        this.end = end;
        this.name = name;
        this.start = start;
        this.timeScores = timeScores;

        float minScore = timeScores[0];
        float maxScore = timeScores[0];
        int minIdx = 0;
        int maxIdx = 0;
        for (int i = 0; i < timeScores.length; i++) {
            if (timeScores[i] < minScore) {
                minIdx = i;
                minScore = timeScores[i];
            }
            if (timeScores[i] > maxScore) {
                maxIdx = i;
                maxScore = timeScores[i];
            }
        }

        foldChange = (maxScore + 1) / (minScore + 1);
        dynamicScore =  (maxScore < 30) ?  0 : (float) (Math.log(foldChange) / Globals.log2);
        if(minIdx > maxIdx) {
            dynamicScore = -dynamicScore;
        }
    }

    public String getChr() {
        return chr;
    }

    @Override
    public String getContig() {
        return chr;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public String getName() {
        return name;
    }

    public float getCombinedScore() {
        return combinedScore;
    }

    public float[] getTimeScores() {
        return timeScores;
    }

    public boolean isDynamic() {
        return dynamic;
    }

    public float getDynamicScore() {
        return dynamicScore;
    }

    // Locus score interface

    public void setStart(int start) {
        this.start = start;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public float getScore() {
        return combinedScore;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public LocusScore copy() {
        return this;
    }

    String valueString;

    public String getValueString(double position, WindowFunction windowFunction) {
        if (valueString == null) {
            StringBuffer buf = new StringBuffer();
            buf.append("Combined Score: " + getScore());
            buf.append("<br>Fold change: " + foldChange);
            buf.append("<br>--------------------");
            for (int i = 0; i < timeScores.length; i++) {
                buf.append("<br>Time point " + (i + 1) + ": " + timeScores[i]);
            }

            valueString = buf.toString();
        }
        return valueString;
    }

    public float getFoldChange() {
        return foldChange;
    }
}
