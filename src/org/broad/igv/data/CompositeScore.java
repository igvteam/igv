/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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
package org.broad.igv.data;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.WindowFunction;

import java.util.List;

/**
 * @author jrobinso
 */
public class CompositeScore implements LocusScore {

    List<LocusScore> scores;
    float score;
    private int start;
    private int end;

    /**
     * Constructs ...
     * mean("Mean"),
     * median("Median"),
     * min("Minimum"),
     * max("Maximum"),
     * percentile10("10th Percentile"),
     * percentile90("90th Percentile"),
     * percentile98("98th Percentile"),
     * stddev("Standard Deviation"),
     * count("Count"),
     * density("Density");
     *
     * @param scores
     * @param windowFunction
     */
    public CompositeScore(List<LocusScore> scores, WindowFunction windowFunction) {
        this.scores = scores;
        start = Integer.MAX_VALUE;
        end = -Integer.MAX_VALUE;
        for (LocusScore score : scores) {
            start = Math.min(start, score.getStart());
            end = Math.max(end, score.getEnd());
        }

        float[] data = new float[scores.size()];
        for (int i = 0; i < scores.size(); i++) {
            data[i] = scores.get(i).getScore();
        }
        score = ProcessingUtils.computeStat(data, windowFunction);
    }

    /**
     * Constructs ...
     *
     * @param sc
     */
    public CompositeScore(CompositeScore sc) {
        this.scores = sc.scores;
        this.start = sc.start;
        this.end = sc.end;
        this.score = sc.score;
    }


    public CompositeScore copy() {
        return new CompositeScore(this);
    }


    public float getScore() {

        // Todo - implementation
        return score;
    }


    public void setConfidence(float confidence) {

        // ignore
    }


    public float getConfidence() {
        return 1.0f;
    }


    public String getValueString(double position, WindowFunction windowFunction) {
        if (scores == null) {
            return "";
        }
        StringBuffer buf = new StringBuffer();
        buf.append("Composite value = " + score + " (" + windowFunction + ")<br>");
        buf.append("-------------------------------<br>");
        buf.append("Components: <br>");
        for (LocusScore s : scores) {
            buf.append(s.getValueString(position, windowFunction));
            buf.append("<br>");

        }
        return buf.toString();
    }

    public String getChr() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }


    public int getStart() {
        return start;
    }


    public void setStart(int start) {
        this.start = start;
    }


    public int getEnd() {
        return end;
    }

    /**
     * Method description
     *
     * @param end
     */
    public void setEnd(int end) {
        this.end = end;
    }

    public int getExtendedStart() {
        return getStart();
    }

    public int getExtendedEnd() {
        return getEnd();
    }
}
