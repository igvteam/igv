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
