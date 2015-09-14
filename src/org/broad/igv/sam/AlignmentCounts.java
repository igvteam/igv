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

package org.broad.igv.sam;

import htsjdk.tribble.Feature;

/**
 * @author Jim Robinson
 * @date 11/22/11
 */
public interface AlignmentCounts extends Feature {

    void incCounts(Alignment alignment);

    int getTotalCount(int pos);

    int getTotalQuality(int pos);

    int getCount(int pos, byte b);

    int getNegCount(int pos, byte b);

    int getPosCount(int pos, byte b);

    int getDelCount(int pos);

    int getInsCount(int pos);

    int getQuality(int pos, byte b);

    int getNumberOfPoints();

    int getMaxCount(int origin, int end);

    String getValueStringAt(int pos);

    boolean isMismatch(int pos, byte ref, String chr, float snpThreshold);

    BisulfiteCounts getBisulfiteCounts();

    void finish();

}
