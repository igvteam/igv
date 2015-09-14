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
