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

package org.broad.igv.event;

/**
 * Events corresponding to a change in viewed area (chromosome, position, and/or zoom).
 * <p>
 * {@code Cause} derived events
 * should cause the data model (e.g. ReferenceFrame) to change their position, this will
 * typically be dispatched by a UI Component.
 * <p>
 * {@code Result} derived events should be dispatched after objects have changed
 * their position, typically to tell UI elements they should repaint
 * User: jacob
 * Date: 2013-Jan-30
 */
public class ViewChange {

    public enum Type {Change, ChromosomeChange, LocusChange}

    boolean recordHistory = false;
    public Type type;
    public String chrName;
    public double start;
    public double end;

    private ViewChange(Type type) {
        this.type = type;
    }

    private ViewChange(Type type, String chrName) {
        this.type = type;
        this.chrName = chrName;
    }

    private ViewChange(Type type, String chrName, double start, double end) {
        this.type = type;
        this.chrName = chrName;
        this.start = start;
        this.end = end;
    }

    public boolean recordHistory() {
        return this.recordHistory;
    }

    public void setRecordHistory(boolean recordHistory) {
        this.recordHistory = recordHistory;
    }


    public static ViewChange Result() {
        return new ViewChange(Type.Change);
    }

    public static ViewChange ChromosomeChangeResult(String chrName) {
        return new ViewChange(Type.ChromosomeChange, chrName);
    }

    public static ViewChange LocusChangeResult(String chrName, double start, double end) {
        return new ViewChange(Type.LocusChange, chrName, start, end);
    }

}
