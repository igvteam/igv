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
 * Events which either cause or are the result of changes in the zoom level
 * User: jacob
 * Date: 2013-Jan-30
 */
public class ZoomChange {


    /**
     * Event indicating that the zoom should change. This event
     * will generally be sent by UI components which want to change
     * the zoom
     */
    public static class Cause{
        //public final int oldZoom;
        public final int newZoom;
        //public final Object source;

        public Cause(int newZoom){
            //this.oldZoom = oldZoom;
            this.newZoom = newZoom;
            //this.source = source;
        }
    }

    /**
     * Event dispatched after objects have changed their zoom level,
     * generally so UI components can repaint
     */
    public static class Result{
        public final int currentZoom;

        public Result(int currentZoom){
            this.currentZoom = currentZoom;
        }
    }
}
