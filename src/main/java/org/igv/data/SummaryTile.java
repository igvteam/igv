package org.igv.data;


import org.igv.feature.LocusScore;

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
