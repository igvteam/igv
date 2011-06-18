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

package org.broad.igv.data.seg;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.WindowFunction;

/**
 * @author Enter your name here...
 * @version Enter version here..., 09/01/09
 */
public class Segment implements LocusScore {

    private int extendedStart = -1;
    private int extendedEnd = -1;
    private int start;
    private int end;
    private float score;
    private String description;


    public Segment(int start, int end, float score) {
        this.start = start;
        this.end = end;
        if (extendedStart < 0) {
            extendedStart = start;
        }
        if (extendedEnd < 0) {
            extendedEnd = end;
        }
        this.score = score;
    }



    public Segment(int start, int origStart, int end, int origEnd, float value, String description) {
        this.start = start;
        this.end = end;
        this.extendedStart = origStart;
        this.extendedEnd = origEnd;
        this.score = value;
        this.description = description;
    }

    public Segment copy() {
        return new Segment(start, extendedStart, end, extendedEnd, score, description);
    }

    public String getChr() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public float getScore() {
        return score;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public void setEnd(int end) {
        this.end = end;
    }


    public String getValueString(double position, WindowFunction ignored) {
        String valueString = "Value: " + getScore();
        if (description != null) {
            valueString += description;
        }
        return valueString;
    }

    public String getDescription() {
        return description;
    }
}
