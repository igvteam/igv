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
