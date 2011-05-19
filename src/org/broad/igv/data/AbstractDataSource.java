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
package org.broad.igv.data;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.tdf.Bin;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.util.LRUCache;
import org.broad.igv.util.collections.DoubleArrayList;
import org.broad.igv.util.collections.FloatArrayList;

import java.util.*;

/**
 * @author jrobinso
 */
public abstract class AbstractDataSource implements DataSource {

    private static Logger log = Logger.getLogger(AbstractDataSource.class);

    // DataManager dataManager;
    boolean cacheSummaryTiles = true;
    WindowFunction windowFunction = WindowFunction.mean;
    LRUCache<String, SummaryTile> summaryTileCache = new LRUCache(this, 20);
    Genome genome;

    public AbstractDataSource(Genome genome) {
        this.genome = genome;
    }

    // abstract DataManager getDataManager();

    /**
     * Return the number of precomputed zoom levels for the given chromosome.
     *
     * @param chr
     * @return
     */
    abstract protected int getNumZoomLevels(String chr);

    // abstract protected TrackType getTrackType();

    /**
     * Return "raw" (i.e. not summarized) data for the specified interval.
     *
     * @param chr
     * @param startLocation
     * @param endLocation
     * @return
     */
    abstract protected DataTile getRawData(String chr, int startLocation, int endLocation);

    /**
     * Return the precomputed summary tiles for the given locus and zoom level.  If
     * there are non return null.
     *
     * @param chr
     * @param startLocation
     * @param endLocation
     * @param zoom
     * @return
     */
    protected List<SummaryTile> getPrecomputedSummaryTiles(String chr, int startLocation, int endLocation, int zoom) {
        return null;
    }


    public int getChrLength(String chr) {
        if (chr.equals(Globals.CHR_ALL)) {
            return (int) (genome.getLength() / 1000);
        } else {
            Chromosome c = genome.getChromosome(chr);
            return c == null ? 0 : c.getLength();
        }
    }

    /**
     * Refresh the underlying data. Default implementation does nothing, subclasses
     * can override
     *
     * @param timestamp
     */
    public void refreshData(long timestamp) {

        // ignore --
    }

    /**
     * Return the longest feature in the dataset for the given chromosome.  This
     * is needed when computing summary data for a region.
     * <p/>
     * TODO - This default implementaiton is crude and should be overriden by subclasses.
     *
     * @param chr
     * @return
     */
    public abstract int getLongestFeature(String chr);

    //{
    //
    //    if (getTrackType() == TrackType.GENE_EXPRESSION) {
    //        String genomeId = IGV.getInstance().getGenomeManager().getGenomeId();
    //        GeneManager gm = GeneManager.getGeneManager(genomeId);
    //        return (gm == null) ? 1000000 : gm.getLongestGeneLength(chr);
    //    } else {
    //        return 1000;
    //    }
    //}


    public List<LocusScore> getSummaryScoresForRange(String chr, int startLocation, int endLocation, int zoom) {


        List<SummaryTile> tiles = getSummaryTilesForRange(chr, startLocation, endLocation, zoom);

        List<LocusScore> summaryScores = new ArrayList(tiles.size() * 700);

        for (SummaryTile tile : tiles) {
            summaryScores.addAll(tile.getScores());

        }
        //FeatureUtils.sortFeatureList(summaryScores);
        return summaryScores;

    }

    private List<SummaryTile> getSummaryTilesForRange(String chr, int startLocation, int endLocation, int zReq) {
        int startTile;
        int endTile;
        if (zReq < getNumZoomLevels(chr)) {
            return getPrecomputedSummaryTiles(chr, startLocation, endLocation, zReq);
        } else {
            int chrLength = getChrLength(chr);
            if (chrLength == 0) {
                return Collections.emptyList();
            }
            endLocation = Math.min(endLocation, chrLength);

            // By definition there are 2^z tiles per chromosome, and 700 bins per tile, where z is the zoom level.
            //int maxZoom = (int) (Math.log(chrLength/700) / Globals.log2) + 1;
            //int z = Math.min(zReq, maxZoom);
            int z = zReq;
            int nTiles = (int) Math.pow(2, z);
            //double binSize = Math.max(1, (((double) chrLength) / nTiles) / 700);


            int adjustedStart = Math.max(0, startLocation);
            int adjustedEnd = Math.min(chrLength, endLocation);


            if (cacheSummaryTiles && !FrameManager.isGeneListMode()) {
                double tileWidth = ((double) chrLength) / nTiles;

                startTile = (int) (adjustedStart / tileWidth);
                endTile = (int) (Math.min(chrLength, adjustedEnd) / tileWidth) + 1;
                List<SummaryTile> tiles = new ArrayList(nTiles);
                for (int t = startTile; t <= endTile; t++) {
                    int tileStart = (int) (t * tileWidth);
                    int tileEnd = Math.min(chrLength, (int) ((t + 1) * tileWidth));

                    String key = chr + "_" + z + "_" + t + getWindowFunction();
                    SummaryTile summaryTile = summaryTileCache.get(key);
                    if (summaryTile == null) {

                        summaryTile = computeSummaryTile(chr, t, tileStart, tileEnd, 700);

                        if (cacheSummaryTiles && !FrameManager.isGeneListMode()) {
                            synchronized (summaryTileCache) {
                                summaryTileCache.put(key, summaryTile);
                            }
                        }
                    }


                    if (summaryTile != null) {
                        tiles.add(summaryTile);
                    }
                }
                return tiles;
            } else {
                SummaryTile summaryTile = computeSummaryTile(chr, 0, startLocation, endLocation, 700);
                return Arrays.asList(summaryTile);
            }


        }
    }

    /**
     * Note:  Package scope used so this method can be unit tested
     *
     * @param chr
     * @param tileNumber
     * @param startLocation
     * @param endLocation
     * @return
     */

    SummaryTile computeSummaryTile(String chr, int tileNumber, int startLocation, int endLocation, int nBins) {


        // TODO -- we should use an index here
        int longestGene = getLongestFeature(chr);

        int adjustedStart = Math.max(startLocation - longestGene, 0);
        DataTile rawTile = getRawData(chr, adjustedStart, endLocation);
        SummaryTile tile = null;

        if ((rawTile != null) && !rawTile.isEmpty() && nBins > 0) {
            int[] starts = rawTile.getStartLocations();
            int[] ends = rawTile.getEndLocations();
            float[] values = rawTile.getValues();
            String[] features = rawTile.getFeatureNames();

            tile = new SummaryTile(tileNumber, startLocation);

            if (windowFunction == WindowFunction.none) {

                for (int i = 0; i < starts.length; i++) {
                    int s = starts[i];
                    int e = ends == null ? s + 1 : Math.max(s + 1, ends[i]);

                    if (e < startLocation) {
                        continue;
                    } else if (s > endLocation) {
                        break;
                    }

                    String probeName = features == null ? null : features[i];
                    float v = values[i];

                    BasicScore score = new BasicScore(s, e, v);
                    tile.addScore(score);

                }


            } else {


                /*
                // Physical overlap
                // TODO -- most datasets do not have data features that physically overlap.  Determine during parsing (not here)
                List<LocusScore> scores = new ArrayList();
                List<PendingScore> pending = new LinkedList();
                for (int i = 0; i < starts.length; i++) {
                    int s = starts[i];
                    int e = ends == null ? s + 1 : Math.max(s + 1, ends[i]);
                    String probeName = features == null ? null : features[i];
                    float v = values[i];

                    if (e < startLocation) {
                        continue;
                    } else if (s > endLocation) {
                        break;
                    }
                    
                    if (pending.isEmpty()) {
                        pending.add(new PendingScore(s, e, v));

                    } else {

                        Iterator<PendingScore> iter = pending.iterator();
                        PendingScore new1 = null;
                        PendingScore new2 = null;
                        while (iter.hasNext()) {
                            PendingScore ps = iter.next();
                            if (ps.end <= s) {
                                scores.add(createScore(ps));
                                iter.remove();

                            }
                            // Must at least overlap, overlap or contained?
                            else if (e >= ps.end) {
                                // overlap.  Cut pending score into 3
                                scores.add(createScore(ps));
                                iter.remove();

                                new1 = (new PendingScore(s, ps.end, ps.scores, v));
                                new2 = (new PendingScore(ps.end, e, v));


                            } else {
                                // contained
                                scores.add(createScore(ps));
                                iter.remove();

                                new1 = (new PendingScore(s, e, ps.scores, v));
                                new2 = (new PendingScore(e, ps.end, ps.scores));

                            }
                        }
                        if (new1 != null) {
                            pending.add(new1);
                        }
                        if (new2 != null) {
                            pending.add(new2);
                        }
                    }

                }

                // Empty "pending" collection
                for (PendingScore ps : pending) {
                    if (ps.scores.size() == 1) {
                        scores.add(new BasicScore(ps.start, ps.end, ps.scores.get(0)));
                    } else {
                        scores.add(new CompositeScore(ps.start, ps.end, ps.scores.toArray(), windowFunction));
                    }

                }

                */
                FloatArrayList accumulatedValues = new FloatArrayList();
                List<String> probes = new ArrayList();
                int accumulatedStart = -1;
                int accumulatedEnd = -1;
                int lastEndBin = 0;

                for (int i = 0; i < starts.length; i++) {
                    int s = starts[i];
                    int e = ends == null ? s + 1 : Math.max(s + 1, ends[i]);
                    String probeName = features == null ? null : features[i];
                    float v = values[i];

                    if (e < startLocation) {
                        continue;
                    } else if (s > endLocation) {
                        break;
                    }


                    double scale = (double) (endLocation - startLocation) / nBins;


                    // if (score.getEnd() < startLocation) continue;
                    // if (score.getStart() > endLocation) break;

                    int startBin = Math.max(0, (int) ((s - startLocation) / scale));
                    int endBin = Math.max(0, (int) ((e - startLocation) / scale));
                    if ((endBin - startBin) > 1) {
                        // This feature covers more than one pixel.  Record previous bin, if any and this one as well.
                        if (!accumulatedValues.isEmpty()) {
                            LocusScore ls = accumulatedValues.size() == 1 ?
                                    new NamedScore(accumulatedStart, accumulatedEnd, accumulatedValues.get(0), probes.get(0)) :
                                    new CompositeScore(accumulatedStart, accumulatedEnd, accumulatedValues.toArray(),
                                            probes.toArray(new String[0]), windowFunction);
                            tile.addScore(ls);
                            accumulatedValues.clear();
                            probes.clear();
                            accumulatedStart = -1;
                            accumulatedEnd = -1;
                        }

                        tile.addScore(new NamedScore(s, e, v, probeName));


                    } else { //endBin == startBin

                        if (endBin == lastEndBin) {
                            // Add to previous bin
                            accumulatedValues.add(v);
                            probes.add(probeName);
                            if (accumulatedStart < 0) accumulatedStart = s;
                            accumulatedEnd = e;
                        } else {
                            // Halt previous "bin" and start a new one
                            if (!accumulatedValues.isEmpty()) {
                                LocusScore ls = accumulatedValues.size() == 1 ?
                                        new NamedScore(accumulatedStart, accumulatedEnd, accumulatedValues.get(0), probes.get(0)) :
                                        new CompositeScore(accumulatedStart, accumulatedEnd, accumulatedValues.toArray(),
                                                probes.toArray(new String[0]), windowFunction);
                                tile.addScore(ls);
                            }
                            accumulatedValues.clear();
                            probes.clear();

                            accumulatedValues.add(v);
                            probes.add(probeName);
                            accumulatedStart = s;
                            accumulatedEnd = e;
                        }

                    }
                    lastEndBin = endBin;
                }
                if (!accumulatedValues.isEmpty()) {
                    LocusScore ls = accumulatedValues.size() == 1 ?
                            new NamedScore(accumulatedStart, accumulatedEnd, accumulatedValues.get(0), probes.get(0)) :
                            new CompositeScore(accumulatedStart, accumulatedEnd, accumulatedValues.toArray(),
                                    probes.toArray(new String[0]), windowFunction);
                    tile.addScore(ls);
                }

            }
        }

        return tile;
    }


    static class PendingScore {
        int start;
        int end;
        FloatArrayList scores = new FloatArrayList();

        PendingScore(int start, int end, FloatArrayList scores, float newVal) {
            this.end = end;
            this.scores.addAll(scores);
            this.scores.add(newVal);
            this.start = start;
        }

        PendingScore(int start, int end, FloatArrayList scores) {
            this.end = end;
            this.scores.addAll(scores);
            this.start = start;
        }

        PendingScore(int start, int end, float newVal) {
            this.end = end;
            this.scores.add(newVal);
            this.start = start;
        }
    }


    /**
     * Return true if the data has been log normalized.
     *
     * @return
     */
    public boolean isLogNormalized() {
        return true;
    }


    public void setWindowFunction(WindowFunction statType) {
        this.windowFunction = statType;
        this.summaryTileCache.clear();
    }


    public WindowFunction getWindowFunction() {
        return windowFunction;
    }

    // TODO -- get window functions dynamically from data
    static List<WindowFunction> wfs = new ArrayList();


    static {
        wfs.add(WindowFunction.min);
        wfs.add(WindowFunction.percentile10);
        wfs.add(WindowFunction.median);
        wfs.add(WindowFunction.mean);
        wfs.add(WindowFunction.percentile90);
        wfs.add(WindowFunction.max);

    }


    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return wfs;
    }


    /*
  SummaryTile computeSummaryTile(String chr, int tileNumber, int startLocation, int endLocation, float binSize) {


      // TODO -- we should use an index here
      int longestGene = getLongestFeature(chr);

      int adjustedStart = Math.max(startLocation - longestGene, 0);
      DataTile rawTile = getRawData(chr, adjustedStart, endLocation);
      SummaryTile tile = null;

      int nBins = (int) ((endLocation - startLocation) / binSize + 1);

      if ((rawTile != null) && !rawTile.isEmpty() && nBins > 0) {
          int[] starts = rawTile.getStartLocations();
          int[] ends = rawTile.getEndLocations();
          float[] values = rawTile.getValues();
          String[] features = rawTile.getFeatureNames();

          tile = new SummaryTile(tileNumber, startLocation);

          if (windowFunction == WindowFunction.none) {

              for (int i = 0; i < starts.length; i++) {
                  int s = starts[i];
                  int e = ends == null ? s + 1 : Math.max(s + 1, ends[i]);

                  if (e < startLocation) {
                      continue;
                  } else if (s > endLocation) {
                      break;
                  }

                  String probeName = features == null ? null : features[i];
                  float v = values[i];

                  BasicScore score = new BasicScore(s, e, v);
                  tile.addScore(score);

              }


          } else {
              List<Bin> bins = new ArrayList();

              int startIdx = 0;
              int endIdx = starts.length;
              for (int i = 0; i < starts.length; i++) {
                  int s = starts[i];
                  int e = ends == null ? s + 1 : Math.max(s + 1, ends[i]);
                  if (e < startLocation) {
                      startIdx = i;
                      continue;
                  } else if (s > endLocation) {
                      endIdx = i;
                      break;
                  }
              }


              for (int i = startIdx; i < endIdx; i++) {
                  int s = starts[i];
                  int e = ends == null ? s + 1 : Math.max(s + 1, ends[i]);
                  float v = values[i];
                  String probeName = features == null ? null : features[i];

                  int start = binPosition(binSize, s);
                  int end = binPosition(binSize, e);

                  // Previous bins are already sliced into non-overlapping section.  The current bin can (1) contribute
                  // to previous bins, (2) slice a previous bin if end intersects it, or (3) start a new bin
                  // if start is >= the end of the last previous bin

                  int lastEnd = 0;
                  int binIdx = 0;
                  for (Bin b : bins) {
                      int bEnd = b.getEnd();
                      lastEnd = bEnd;
                      if (b.getStart() >= end) {
                          // This should never happen?
                          break;
                      }
                      if (b.getEnd() <= start) {
                          continue;
                      }

                      if (b.getStart() < start) {
                          //    [     ]            <= b
                          //        [     ]
                          // Shift end of bin
                          b.setEnd(start);

                          // Create the merged bin, representing the overlapped region
                          Bin mergedBin = new Bin(start, bEnd, probeName, v, windowFunction);
                          mergedBin.mergeValues(b);
                          bins.add(mergedBin);

                          // Create a new bin for the non-overlapped part
                          Bin newBin = new Bin(bEnd, end, probeName, v, windowFunction);
                          bins.add(newBin);

                      } else {
                          //     [     ]          <= b
                          //       [ ]
                          // Shift end of bin
                          b.setEnd(start);

                          // Create the merged bin, representing the overlapped region
                          Bin mergedBin = new Bin(start, bEnd, probeName, v, windowFunction);
                          mergedBin.mergeValues(b);
                          bins.add(mergedBin);

                          // Create a new bin for the non-overlapped part
                          Bin newBin = new Bin(bEnd, end, windowFunction);
                          newBin.mergeValues(b);
                          bins.add(newBin);
                      }
                      binIdx++;
                  }

                  if (start >= lastEnd) {
                      bins.add(new Bin(start, end, probeName, v, windowFunction));
                  }

              }

              tile.addAllScores(bins);
          }
      }

      return tile;
  }

  private int binPosition(float binSize, int s) {
      return (int) (Math.round(s / binSize) * binSize);
  }  */

}
