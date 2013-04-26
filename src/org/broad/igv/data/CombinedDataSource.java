/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
 * Data source which combines two other DataSources
 *
 * TODO Multiple DataSources. There is no 2
 * @author jrobinso, jacob
 */
public class CombinedDataSource implements DataSource {

    enum Operation{
        ADD("+"),
        SUBTRACT("-");

        private String stringRep;

        private Operation(String stringrep){
            this.stringRep = stringrep;
        }

    }

    DataSource source1;
    DataSource source2;

    Operation operation = Operation.ADD;

    CombinedDataSource(DataSource source1, DataSource source2, Operation operation){
        this.source1 = source1;
        this.source2 = source2;
        this.operation = operation;
    }

    public List<LocusScore> getSummaryScoresForRange(String chr, int startLocation, int endLocation, int zoom){

        List<LocusScore> outerScores = this.source1.getSummaryScoresForRange(chr, startLocation, endLocation, zoom);
        List<LocusScore> innerScores = this.source2.getSummaryScoresForRange(chr, startLocation, endLocation, zoom);

        List<LocusScore> combinedScores = new ArrayList<LocusScore>(outerScores.size() + innerScores.size());

        if(combinedScores.size() == 0) return combinedScores;

        //TODO We assume that having no data from one source is the identity operation, that may not be true
        if(innerScores.size() == 0) return outerScores;
        if(outerScores.size() == 0) return innerScores;

        /**
         * Pretty simple algorithm.  Maybe not super-efficient
         * TODO Check efficiency
         * Assume the following tiles (no breaks)
         * outerScores: |----|--------|--------------|-------|-----------|
         * innerScores:           |-----|-----|--------|-----------|-----------|---------|
         * combined:    |----|----|---|-|-----|------|-|-----|-----|-----|-----|---------|
         *
         * We loop through outerScores, and on each loop through innerScores
         * The overlap region between each is identified, and the scores combined
         */


        //We must have the outerScores variable to start not later than innerScores
        if(outerScores.get(0).getStart() > innerScores.get(0).getStart()){
            List<LocusScore> tmp = innerScores;
            innerScores = outerScores;
            outerScores = tmp;
        }

        int firstInnerStart = innerScores.get(0).getStart();
        int firstInnerEnd = innerScores.get(0).getEnd();

        int maxCoord = 0, highestInnerIdx = -1;

        for(LocusScore outerScore: outerScores){

            int outerStart = outerScore.getStart();
            int outerEnd = outerScore.getEnd();
            highestInnerIdx = -1;

            //Add in regions where outerScores has data but innerScores doesn't
            if(firstInnerStart > outerEnd){
                combinedScores.add(new BasicScore(outerStart, outerEnd, combineScores(outerScore, null)));
                continue;
            }

            for(LocusScore innerScore: innerScores){
                //Past the overlapping region
                if(innerScore.getStart() >= outerEnd) break;

                highestInnerIdx++;

                //Have not yet reached overlapping region
                if(innerScore.getEnd() < outerStart) continue;


                int combinedEnd = Math.min(outerEnd, innerScore.getEnd());
                float combinedVal = combineScores(outerScore, innerScore);
                BasicScore newScore = new BasicScore(outerStart, combinedEnd, combinedVal);

                combinedScores.add(newScore);
                maxCoord = combinedEnd;

            }
        }

        //Get the remaining innerScores
        LocusScore innerTail = innerScores.get(highestInnerIdx);
        int combinedStart = Math.min(maxCoord, innerTail.getEnd());
        BasicScore newTail = new BasicScore(combinedStart, innerTail.getEnd(), innerTail.getScore());
        combinedScores.add(newTail);
        for(LocusScore innerScore: innerScores.subList(highestInnerIdx+1, innerScores.size())){
            float newVal = combineScores(null, innerScore);
            combinedScores.add(new BasicScore(innerScore.getStart(), innerScore.getEnd(), newVal));
        }

        return combinedScores;
    }

    /**
    public List<LocusScore> getSummaryScoresForRange(String chr, int startLocation, int endLocation, int zoom) {

        List<LocusScore> scores1 = getSummaryScoresForRange(chr, startLocation, endLocation, zoom);
        List<LocusScore> scores2 = getSummaryScoresForRange(chr, startLocation, endLocation, zoom);

        List<LocusScore> pendingScores = new ArrayList<LocusScore>();
        List<LocusScore> scoresToClose = new ArrayList<LocusScore>();
        List<LocusScore> mergedScores = new ArrayList<LocusScore>(scores1.size() + scores2.size());

        Iterator<LocusScore> iter = new MergedIterator(scores1.iterator(), scores2.iterator());
        while (iter.hasNext()) {

            LocusScore score = iter.next();
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
     **/

    /**
     * Combine the scores using this sources operation. Either can be null,
     * in which case it is given the coordinates of the other and a score of 0.
     * Both inputs cannot be null
     * @param score1
     * @param score2
     * @return
     */
    private float combineScores(LocusScore score1, LocusScore score2) {
        if(score1 == null && score2 == null) throw new IllegalArgumentException("Both inputs cannot be null");
        if(score1 == null){
            score1 = new BasicScore(score2.getStart(), score2.getEnd(), 0.0f);
        }else if(score2 == null){
            score2 = new BasicScore(score1.getStart(), score1.getEnd(), 0.0f);
        }

        switch(operation){
            case ADD:
                return score1.getScore() + score2.getScore();
            case SUBTRACT:
                return score1.getScore() - score2.getScore();
            default:
                throw new IllegalStateException("Operation not recognized: " + operation);
        }
    }

    /**
     * Iterator which returns combined iterators sorted by start position.
     * i1 comes first in case of ties
     */
    static class MergedIterator implements Iterator<LocusScore> {

        Iterator i1;
        Iterator i2;
        LocusScore next1;
        LocusScore next2;

        MergedIterator(Iterator<LocusScore> i1, Iterator<LocusScore> i2) {
            this.i1 = i1;
            this.i2 = i2;

            if (i1.hasNext()) {
                next1 = i1.next();
            }

            if (i2.hasNext()) {
                next2 = i2.next();
            }
        }

        public boolean hasNext() {
            return next1 != null || next2 != null;
        }

        public LocusScore next() {
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
