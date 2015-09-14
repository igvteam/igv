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

import org.broad.igv.feature.LocusScore;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Represents a chunk of data over a fixed window span (default == 1)
 *
 * @version 1.0
 * @created 13-Nov-2007 13:11:15
 */
public class SummaryTile2D {

    /**
     * The tileNumber.  Used to support unit tests, and  general  debugging
     */
    private int tileNumber;
    private int startLocation;

    Map<Integer, List<LocusScore>> summaryScores;

    public SummaryTile2D(int tileNumber, int startLocation) {
        this.tileNumber = tileNumber;
        this.startLocation = startLocation;
        summaryScores = new HashMap();
    }

    public void addScore(int trackNumber, LocusScore score) {
        List<LocusScore> scoreList = summaryScores.get(trackNumber);
        if (scoreList == null) {
            scoreList = new ArrayList();
            summaryScores.put(trackNumber, scoreList);
        }
        scoreList.add(score);
    }

    public List<LocusScore> getScores(int trackNumber) {
        return summaryScores.get(trackNumber);
    }

    public Map<Integer, List<LocusScore>> getSummaryScores() {
        return summaryScores;
    }

    public int getSize() {
        return summaryScores.size();
    }

    public boolean isEmpty() {
        return summaryScores.isEmpty();
    }

    public int getStartLocation() {
        return startLocation;
    }


}
