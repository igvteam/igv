/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

/*
 * To change this template, choose Tools | Templates
 * and openFile the template in the editor.
 */
package org.broad.igv.data;

import org.broad.igv.Globals;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * @author jrobinso
 */
public class HDFDataSource implements DataSource {

    private String name;
    private int trackNumber;  // Order in HDF5File

    //private int version;
    HDFDataManager dataManager;

    public HDFDataSource(HDFDataManager dataManager, String name, int trackNumber) {

        this.dataManager = dataManager;
        this.name = name;
        this.trackNumber = trackNumber;

    }

    public boolean isCopyNumber() {
        TrackType tt = dataManager.getTrackType();
        return tt == TrackType.COPY_NUMBER || tt == TrackType.ALLELE_SPECIFIC_COPY_NUMBER;
    }

    String getName() {
        return name;
    }

    public TrackType getTrackType() {
        return dataManager.getTrackType();
    }

    public void setWindowFunction(WindowFunction statType) {
        dataManager.setWindowFunction(statType);
    }

    public boolean isLogNormalized() {
        return dataManager.isNormalized();
    }

    public int getChrLength(String chr) {
        return dataManager.getChrLength(chr);
    }

    public double getDataMax() {
        return dataManager.getDataMax(trackNumber, Globals.CHR_ALL);
    }

    public double getDataMin() {
        return dataManager.getDataMin(trackNumber, Globals.CHR_ALL);
    }

    /**
     * Refresh the data
     */
    public void refresh() {
        // TODO -- implementation.  Clear caches, etc.
    }

    /**
     * @param trackNumber
     * @param chr
     * @param startLocation
     * @param endLocation
     * @param zoom
     * @param windowFunction
     * @return
     */
    // TDO -- without caching the lastTileList the app hangs in some rare cases
    //  FIGURE OUT WHY!!!
    private static String lastTileListKey = null;
    private static List<SummaryTile2D> lastTileList = null;

    public List<LocusScore> getSummaryScoresForRange(
            String chr, int startLocation, int endLocation, int zoom) {

        List<SummaryTile2D> tiles = null;
        String key = dataManager.getResourceLocator().toString() + "_" +
                chr + " " + startLocation + "_" + endLocation + "_" +
                zoom + getWindowFunction().toString();
        if (lastTileListKey != null && lastTileListKey.equals(key)) {
            tiles = lastTileList;
        } else {
            tiles = dataManager.getSummaryTilesForRange(chr, startLocation, endLocation, zoom);
            lastTileListKey = key;
            lastTileList = tiles;
        }


        List<LocusScore> summaryScores = new ArrayList(tiles.size() * 700);

        for (SummaryTile2D tile : tiles) {
            if (!tile.isEmpty()) {
                List<LocusScore> scores = tile.getScores(trackNumber);
                if (scores != null) {
                    summaryScores.addAll(tile.getScores(trackNumber));
                }
            }
        }

        return summaryScores;

    }


    /**
     * Method supporting obsolete option to join adjacent copy number scores.
     *
     * @param trackNumber
     * @param tiles
     * @return
     * @deprecated Option no longer supported
     */
    private List<LocusScore> aggregateScores(int trackNumber, List<SummaryTile2D> tiles) {
        List<LocusScore> joinedScores = new ArrayList(tiles.size() * 700);
        LocusScore previousScore = null;

        for (SummaryTile2D tile : tiles) {
            if (tile.getScores(trackNumber) != null) {
                for (LocusScore score : tile.getScores(trackNumber)) {
                    if (!Float.isNaN(score.getScore())) {

                        if (previousScore == null) {
                            previousScore = new SummaryScore(score);
                            joinedScores.add(previousScore);
                        } else {
                            if (score.getScore() == previousScore.getScore()) {
                                // score values are identical, stretch the current score
                                // The only purpose of this is to reconsitute segmented data from
                                // cn files that were created from cbs (segmented) files.  This
                                // really doesn'tileNumber belong in IGV.
                                previousScore.setEnd(score.getEnd());
                                previousScore.setConfidence(1);

                            } else {

                                SummaryScore newScore = new SummaryScore(score);

                                // score value has changed. Adjust end of previous
                                // score (if any), and start of this score to meet 1/2
                                // way
                                int delta = newScore.getStart() - previousScore.getEnd();
                                previousScore.setEnd(previousScore.getEnd() + delta / 2);
                                newScore.setStart(previousScore.getEnd());

                                joinedScores.add(newScore);

                                previousScore = newScore;
                            }

                        }
                    }
                }

            }
        }


        return joinedScores;
    }

    public void refreshData(long timestamp) {
        // ignored for now
    }

    public Collection<WindowFunction> getAvailableWindowFunctions() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public WindowFunction getWindowFunction() {
        return dataManager.getWindowFunction();
    }
}
