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

import org.broad.igv.track.RenderContext;

/**
 * User: jacob
 * Date: 2013-Feb-06
 */
public class DataLoadedEvent {

    public final RenderContext context;

    public DataLoadedEvent(RenderContext context){
        this.context = context;
    }
}
