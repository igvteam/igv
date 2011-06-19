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
