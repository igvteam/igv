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
