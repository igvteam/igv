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
 * and open the template in the editor.
 */
package org.broad.igv.tdf;

import org.broad.igv.track.WindowFunction;

import java.util.List;

/**
 * @author jrobinso
 */
public class TileManager {

    TDFReader reader;

    public TileManager(TDFReader reader) {
        this.reader = reader;
    }

    public List<TDFTile> getTiles(String chr, int start, int end, int zoom) {

        TDFDataset ds = reader.getDataset(chr, zoom, WindowFunction.mean);
        if (ds != null) {
            return ds.getTiles(start, end);
        }
        return null;
    }

    /*
    private TDFTile computeTile(String chr, int start, int end, int zoom) {
        String dsName = "/" + chr + "/raw";

        TDFDataset ds = reader.getDataset(dsName);
        if (ds != null) {

            double binSize = tileWidth / 700;
            int binSizeInt = (int) binSize;

            List<TDFTile> tiles = ds.getTiles(startLocation, endLocation);
            if (tiles.size() > 0) {

                for (TDFTile tile : tiles) {
                    // Tile of raw data
                    if (tile != null && tile.getSize() > 0) {

                        int lastBin = -1;
                        int nPts = 0;
                        float[] sum = new float[nTracks];
                        for (int i = 0; i < tile.getSize(); i++) {
                            int s = tile.getStartPosition(i);
                            int e = tile.getEndPosition(i);
                            int bin = (int) (s / binSize);

                            for (int t = 0; i < nTracks; t++) {

                                // TODO -- that overlap and span multiple bins
                                if ((e - s) >= binSizeInt) {
                                } else {
                                    if (lastBin < 0 || bin == lastBin) {
                                        sum[t] = sum[t] + tile.getValue(t, i);
                                        nPts++;
                                    } else {
                                        // On to a new bin.  Record previous one and start over
                                        // Note we're only doing the mean here,  should use the "window function"
                                        float mean = sum[t] / nPts;
                                        int binStart = (int) (lastBin * binSize);
                                        int binEnd = (int) Math.max(binStart + 1, (lastBin + 1) * binSize);
                                        //scores.add(new BasicScore(chr, binStart, binEnd, mean));

                                        sum[t] = v[t];
                                        nPts = 1;
                                    }
                                }
                                lastBin = bin;
                            }
                        }
                    }
                }
            }
        }
    }
     * */
}
