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
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.data;

import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TrackType;

/**
 * A dataset is an in-memory representation of a numerical dataset.  It is used for non-indexed data formats that are
 * read into memory as opposed to read off disk as needed.   The use of objects is eschewed in favor of simple arrays
 * to minimize memory.  
 *
 *
 * @author jrobinso
 */
public interface Dataset {

    public String getName();

    public TrackType getType();

    public TrackProperties getTrackProperties();
    
    public float getDataMin();

    public float getDataMax();

    public String[] getChromosomes();

    public String[] getTrackNames();

    public int[] getStartLocations(String chr);

    public int[] getEndLocations(String chr);

    public String[] getFeatureNames(String chr);

    public float[] getData(String trackName, String chr);

    public boolean isLogNormalized();

    Integer getLongestFeature(String chr);
}
