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

import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

/**
 * @author jrobinso
 * @date Jan 18, 2011
 */
public class CombinedDataSource implements DataSource {

    DataSource source1;
    DataSource source2;
    String operator = "-";   // + or -


    public List<LocusScore> getSummaryScoresForRange(String chr, int startLocation, int endLocation, int zoom) {

        List<LocusScore> scores1 = getSummaryScoresForRange(chr, startLocation, endLocation, zoom);
        List<LocusScore> scores2 = getSummaryScoresForRange(chr, startLocation, endLocation, zoom);

        List<LocusScore> pendingScores = new ArrayList();
        List<LocusScore> scoresToClose = new ArrayList();
        List<LocusScore> mergedScores = new ArrayList(scores1.size() + scores2.size());

        Iterator<ScoreWrapper> iter = new MergedIterator(scores1.iterator(), scores2.iterator());
        while (iter.hasNext()) {

            ScoreWrapper score = iter.next();
            int start = score.getStart();

            // Loop through the scores that are left cutting them at "score" boundary
            for (LocusScore ps : pendingScores) {

                if (ps.getEnd() <= start) {
                    // Done with ps
                    mergedScores.add(ps);
                    scoresToClose.add(ps);

                } else {
                    BasicScore newScore = new BasicScore(ps.getStart(), start, ps.getScore());
                    mergedScores.add(newScore);


                    // The common chunk, we know ps.end is > score.start
                    int end = Math.min(ps.getEnd(), score.getEnd());
                    float newVal = combineScores(ps, score);
                    newScore = new BasicScore(start, end, newVal);
                    pendingScores.add(newScore);

                    // The "tail" of ps,  if any
                    if (ps.getEnd() > score.getEnd()) {
                        newScore = new BasicScore(score.getEnd(), ps.getEnd(), ps.getScore());
                        pendingScores.add(newScore);
                    }
                }

            }
            pendingScores.removeAll(scoresToClose);

        }


        return null;
    }

    private float combineScores(LocusScore score1, ScoreWrapper score2) {
        return score1.getScore() + score2.getScaledScore();
    }


    static class MergedIterator implements Iterator<ScoreWrapper> {

        Iterator i1;
        Iterator i2;
        ScoreWrapper next1;
        ScoreWrapper next2;

        MergedIterator(Iterator<LocusScore> i1, Iterator<LocusScore> i2) {
            this.i1 = i1;
            this.i2 = i2;

            if (i1.hasNext()) {
                next1 = new ScoreWrapper(i1.next(), 1);
            }

            if (i2.hasNext()) {
                next2 = new ScoreWrapper(i1.next(), -1);   // <= subtraction operation
            }
        }

        public boolean hasNext() {
            return next1 != null || next2 != null;
        }

        public ScoreWrapper next() {
            if (next1 == null) {
                return next2;
            } else if (next2 == null) {
                return next1;
            } else if (next2.getStart() < next1.getStart()) {
                return next2;
            } else {
                return next1;
            }
        }

        public void remove() {
            //ignore
        }
    }

    static class ScoreWrapper {
        int operator; // 1 or -1
        LocusScore score;

        ScoreWrapper(LocusScore score, int operator) {
            this.operator = operator;
            this.score = score;
        }

        int getStart() {
            return score.getStart();
        }

        int getEnd() {
            return score.getEnd();
        }

        float getScaledScore() {
            return operator * score.getScore();
        }
    }


    public double getDataMax() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public double getDataMin() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }


    public TrackType getTrackType() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public void setWindowFunction(WindowFunction statType) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public boolean isLogNormalized() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public void refreshData(long timestamp) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public WindowFunction getWindowFunction() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
