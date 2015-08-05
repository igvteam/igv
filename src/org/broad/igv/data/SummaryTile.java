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


package org.broad.igv.data;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.feature.LocusScore;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Represents a chunk of data over a fixed window span (default == 1)
 *
 * @version 1.0
 * @created 13-Nov-2007 13:11:15
 */
public class SummaryTile {

     private int startLocation;

    List<LocusScore> summaryScores;


    public SummaryTile() {
        summaryScores = new ArrayList(1000);
    }

    public SummaryTile(List<LocusScore> summaryScores) {
        this.summaryScores = summaryScores;
    }

    public void addScore(LocusScore score) {
        summaryScores.add(score);
    }


    public void addAllScores(Collection<? extends LocusScore> scores) {
        summaryScores.addAll(scores);
    }


    public List<LocusScore> getScores() {
        return summaryScores;
    }


    public int getSize() {
        return summaryScores.size();
    }


    public boolean isEmpty() {
        return summaryScores.isEmpty();
    }


}
