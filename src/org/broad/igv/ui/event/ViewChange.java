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

package org.broad.igv.ui.event;

/**
 * Events corresponding to a change in viewed area (chromosome, position, and/or zoom).
 *
 * {@code Cause} derived events
 * should cause the data model (e.g. ReferenceFrame) to change their position, this will
 * typically be dispatched by a UI Component.
 *
 * {@code Result} derived events should be dispatched after objects have changed
 * their position, typically to tell UI elements they should repaint
 * User: jacob
 * Date: 2013-Jan-30
 */
public abstract class ViewChange {
    protected boolean recordHistory = false;

    public boolean recordHistory(){
        return this.recordHistory;
    }

    public void setRecordHistory(boolean recordHistory){
        this.recordHistory = recordHistory;
    }

    public static class Cause extends ViewChange{}

    public static class Result extends ViewChange{}

    /**
     * Event indicating that the zoom should change.
     */
    public static class ZoomCause extends Cause{
        //public final int oldZoom;
        public final int newZoom;
        //public final Object source;

        public ZoomCause(int newZoom){
            //this.oldZoom = oldZoom;
            this.newZoom = newZoom;
            //this.source = source;
        }
    }

    public static class ZoomResult extends Result{}

    public static class ChromosomeChangeCause extends Cause{

        public final Object source;
        public final String chrName;

        /**
         *
         * @param source The object which originated the chromosome change
         * @param chrName
         */
        public ChromosomeChangeCause(Object source, String chrName){
            this.source = source;
            this.chrName = chrName;
        }
    }

    public static class ChromosomeChangeResult extends Result{

        public final Object source;
        public final String chrName;

        public ChromosomeChangeResult(Object source, String chrName){
            this.source = source;
            this.chrName = chrName;
        }
    }

    public static class LocusChangeResult extends Result {

        public final String chrName;
        public final double start;
        public final double end;

        public LocusChangeResult(String chrName, double start, double end) {
            this.chrName = chrName;
            this.start = start;
            this.end = end;
        }
    }

}
