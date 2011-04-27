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

package org.broad.igv.tdf;

import org.broad.igv.util.collections.FloatArrayList;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.tools.ListAccumulator;
import org.broad.igv.track.WindowFunction;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Dec 18, 2009
 * Time: 11:22:54 PM
 * To change this template use File | Settings | File Templates.
 */
public class Bin implements LocusScore {

    private int start;
    private int end;

    Accumulator accumulator;
    private FloatArrayList values;
    private List<String> names;
    WindowFunction windowFunction;
    private static final int maxValues = 5;

    public Bin(int start, int end, WindowFunction windowFunction) {
        this.start = start;
        this.end = end;
        this.values = new FloatArrayList(maxValues);
        this.accumulator = new Accumulator(windowFunction);
        this.windowFunction = windowFunction;
    }

    public Bin(int start, int end, String probeName, float initialValue, WindowFunction windowFunction) {
        this.start = start;
        this.end = end;
        this.values = new FloatArrayList(maxValues);
        this.accumulator = new Accumulator(windowFunction);
        this.windowFunction = windowFunction;
        addValue(probeName, initialValue);
    }

    /**
     * Copy constructor
     *
     * @param otherBin
     */
    public Bin(Bin otherBin) {
        this.start = otherBin.start;
        this.end = otherBin.end;
        this.accumulator = otherBin.accumulator;
        this.names = otherBin.names;
        this.windowFunction = otherBin.windowFunction;
    }

    public boolean isExtension(Bin bin) {
        return (end == bin.start) && getScore() == bin.getScore();
    }


    public void addValue(String name, float value) {
        if (values.size() < maxValues) {
            if (name != null) {
                if (names == null) {
                    names = new ArrayList<String>(maxValues);
                }
                names.add(name);
            }
            values.add(value);
        }
        accumulator.add(value);
    }


    public String getChr() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public int getStart() {
        return start;  //To change body of implemented methods use File | Settings | File Templates.
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

    public float getScore() {
        return accumulator.getValue();
    }

    public LocusScore copy() {
        return new Bin(this);  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String getValueString(double position, WindowFunction windowFunction) {
        StringBuffer sb = new StringBuffer(50);
        sb.append("Value: ");
        sb.append(String.valueOf(getScore()));
        if (values.size() == 1) {
            if (names != null && names.size() > 0) {
                sb.append(" (");
                sb.append(names.get(0));
                sb.append(")");
            }
        } else {
            sb.append("<br> ");
            sb.append(windowFunction.getDisplayName() + " of " +
                    (values.size() == maxValues ? "> " : "") + values.size() + " values:");
            float[] v = values.toArray();
            for (int i = 0; i < v.length; i++) {
                sb.append("<br>   " + v[i]);
                if (names != null && names.size() > i) {
                    sb.append("  (");
                    sb.append(names.get(i));
                    sb.append(")");
                }
            }
            if (v.length == maxValues) {
                sb.append("<br>...");
            }
        }

        return sb.toString();
    }

    public FloatArrayList getValues() {
        return values;
    }

    public List<String> getNames() {
        return names;
    }

    public void mergeValues(Bin b) {
        FloatArrayList vList = b.getValues();
        List<String> probes = b.getNames();
        if(vList.size() != probes.size()) {
            // TODO -- error?
        }
        for(int i=0; i<probes.size(); i++) {
            this.addValue(probes.get(i), vList.get(i));
        }
    }
}
