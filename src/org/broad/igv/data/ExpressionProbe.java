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

/*
 * ExpressionProbe.java
 *
 * Created on October 31, 2007, 8:21 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package org.broad.igv.data;

import org.broad.igv.feature.IGVFeature;

/**
 * @author jrobinso
 */
public class ExpressionProbe implements Comparable {


    private String featureName;
    private String probe;
    private String chr = "";
    private int start = 0;
    private int end = 0;

    public ExpressionProbe(String probe) {
        this.probe = probe;
        this.featureName = probe;
    }

    public ExpressionProbe(String probe, String feature) {
        this.probe = probe;
        this.featureName = feature;
    }

    public ExpressionProbe(IGVFeature feature) {
        this(feature.getName());
        this.chr = feature.getChr();
        this.start = feature.getStart();
        this.end = feature.getEnd();
    }

    public String getFeature() {
        return featureName;
    }

    public String getName() {
        return probe;
    }

    public void setProbe(String probe) {
        this.probe = probe;
    }

    public String getChr() {
        return chr;
    }

    public void setChr(String chr) {
        this.chr = chr;
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

    public void setEnd(int end) {
        this.end = end;
    }

    public String toString() {
        return probe + " " + chr + ":" + start + "-" + end;

    }

    public int compareTo(Object anotherProbe) {
        if (anotherProbe instanceof ExpressionProbe) {
            return getStart() - ((ExpressionProbe) anotherProbe).getStart();
        } else return 0;
    }


}
