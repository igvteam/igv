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

package org.broad.igv.peaks;

import org.broad.igv.Globals;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.WindowFunction;
import org.broad.tribble.Feature;

/**
 * Note:  implementing tribble.Feature will allow us to index these files in the future.
 *
 * @author jrobinso
 * @date Apr 22, 2011
 */
public class Peak implements LocusScore, org.broad.tribble.Feature {

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
