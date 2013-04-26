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

package org.broad.igv.data.seg;

import org.broad.igv.data.BasicScore;
import org.broad.igv.track.WindowFunction;

/**
 * @author jacob
 * @date 09/01/09
 */
public class Segment extends BasicScore {

    private int extendedStart = -1;
    private int extendedEnd = -1;
    private String description;


    public Segment(int start, int end, float score) {
        super(start, end, score);
        if (extendedStart < 0) {
            extendedStart = start;
        }
        if (extendedEnd < 0) {
            extendedEnd = end;
        }
    }

    public Segment(int start, int origStart, int end, int origEnd, float value, String description) {
        super(start, end, value);
        this.extendedStart = origStart;
        this.extendedEnd = origEnd;
        this.description = description;
    }

    public Segment copy() {
        return new Segment(start, extendedStart, end, extendedEnd, score, description);
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
