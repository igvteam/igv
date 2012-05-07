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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.feature;

import org.broad.igv.track.WindowFunction;

/**
 * @author jrobinso
 */
public interface LocusScore extends org.broad.tribble.Feature {

    public void setStart(int start);

    public void setEnd(int end);

    public float getScore();

    /**
     * Return a string to be used for popup text.   The WindowFunction is passed
     * in so it can be used to annotate the value.  The LocusScore object itself
     * does not "know" from what window function it was derived
     *
     * @param position       Zero-based genome position
     * @param windowFunction
     * @return
     */
    public String getValueString(double position, WindowFunction windowFunction);
}
