/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
