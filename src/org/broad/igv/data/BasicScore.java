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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.data;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.WindowFunction;

/**
 * @author jrobinso
 */
public class BasicScore implements LocusScore {

    protected int start;
    protected int end;
    protected float score;

    public BasicScore(int start, int end, float score) {

        this.start = start;
        this.end = end;
        this.score = score;
    }

    public BasicScore(BasicScore bs) {
        this.start = bs.start;
        this.end = bs.end;
        this.score = bs.score;
    }

    public BasicScore copy() {
        return new BasicScore(this);
    }

    /**
     * This method is required by the Tribble interface but not used.  To save space the chromosome is not stored, return null
     * @return
     */
    public String getChr() {
        return null;
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

    @Override
    public String toString() {
        return String.format("BasicScore: %d-%d ; %f", getStart(), getEnd(), getScore());
    }

    public String getValueString(double position, WindowFunction windowFunction) {
        StringBuffer buf = new StringBuffer();
        buf.append(String.format("Value: %g at position %d",  score, (int)position));
        if(windowFunction != null) {
            buf.append("<br>Window function: " + windowFunction);
        }
        return buf.toString();
    }

}
