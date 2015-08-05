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
package org.broad.igv.tools.parsers;

import org.broad.igv.track.TrackType;

/**
 * @author jrobinso
 */
public interface DataConsumer {

    public void setType(String type);

    public void addData(String chr, int start, int end, float[] data, String name);

    public void parsingComplete();

    public void setTrackParameters(TrackType trackType, String trackLine, String[] trackNames);

    void setTrackParameters(TrackType trackType, String trackLine, String[] trackNames, boolean b);

    /**
     * Set a tolerance for "sortedness" of the data.  A start position can be less than
     * the immediately previous start position by this amount.  This is needed for
     * chip-seq processing where the start position of an alignment can be artificially
     * extended post sorting.
     */
    public void setSortTolerance(int tolerance);

    public void setAttribute(String key, String value);

}
