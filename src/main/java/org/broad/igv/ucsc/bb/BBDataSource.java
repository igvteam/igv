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

package org.broad.igv.ucsc.bb;

import htsjdk.samtools.util.Locatable;
import org.broad.igv.Globals;
import org.broad.igv.data.AbstractDataSource;
import org.broad.igv.data.BasicScore;
import org.broad.igv.data.DataSource;
import org.broad.igv.data.DataTile;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;

import java.io.IOException;
import java.util.*;


public class BBDataSource extends AbstractDataSource implements DataSource {

    private final Genome genome;
    Collection<WindowFunction> availableWindowFunctions =
            Arrays.asList(WindowFunction.min, WindowFunction.mean, WindowFunction.max, WindowFunction.none);

    BBFile reader;

    private Map<WindowFunction, List<LocusScore>> wholeGenomeScores;

    private double dataMin = 0;
    private double dataMax = 100;

    List<LocusScore> wgScores = null;

    public BBDataSource(BBFile reader, Genome genome) throws IOException {
        super(genome);
        this.reader = reader;
        this.genome = genome;
        this.wholeGenomeScores = new HashMap<>();
        initMinMax();
    }

    /**
     * Set the "min" and "max" from total summary data.
     */
    private void initMinMax() {

        BBTotalSummary totalSummary = reader.getTotalSummary();
        if (totalSummary == null) {
            dataMin = 0;
            dataMax = 100;
        } else {
            dataMin = totalSummary.minVal;
            if (totalSummary.basesCovered < 100) {
                dataMin = Math.min(0, totalSummary.minVal);
                dataMax = totalSummary.maxVal;
            } else if (totalSummary.minVal < 0) {
                dataMax = Math.min(totalSummary.maxVal, 3 * totalSummary.stddev);
                dataMin = Math.max(totalSummary.minVal, -dataMax);
            } else {
                dataMin = 0;
                dataMax = Math.min(totalSummary.maxVal, totalSummary.mean + 2 * totalSummary.stddev);
            }
        }
    }

    /**
     * Set the "min" and "max" from 1MB resolutiond data.  Read a maximum of 10,000 points for this
     */
//    private void initMinMax2() {
//
//        final int oneMB = 1000000;
//        final BBZoomHeader zoomLevelHeader = zoomLevelForScale(oneMB);
//
//        int nValues = 0;
//        double[] values = new double[10000];
//
//        if (zoomLevelHeader == null) {
//            List<String> chrNames = reader.getChromosomeNames();
//            for (String chr : chrNames) {
//                BigWigIterator iter = reader.getBigWigIterator(chr, 0, chr, Integer.MAX_VALUE, false);
//                while (iter.hasNext()) {
//                    WigItem item = iter.next();
//                    values[nValues++] = item.getWigValue();
//                    if (nValues >= 10000) break;
//                }
//            }
//        } else {
//
//            int z = zoomLevelHeader.getZoomLevel();
//            ZoomLevelIterator zlIter = reader.getZoomLevelIterator(z);
//            if (zlIter.hasNext()) {
//                while (zlIter.hasNext()) {
//                    ZoomDataRecord rec = zlIter.next();
//                    values[nValues++] = (rec.getMeanVal());
//                    if (nValues >= 10000) {
//                        break;
//                    }
//                }
//            }
//        }
//
//        if (nValues > 0) {
//            dataMin = StatUtils.percentile(values, 0, nValues, 10);
//            dataMax = StatUtils.percentile(values, 0, nValues, 90);
//        } else {
//            dataMin = 0;
//            dataMax = 100;
//        }
//    }
    public double getDataMax() {
        return dataMax;
    }

    public double getDataMin() {
        return dataMin;
    }

    @Override
    protected DataTile getRawData(String chr, int start, int end) {

        try {
            long rTreeOffset = reader.getHeader().fullIndexOffset;
            List<byte[]> chunks = this.reader.getLeafChunks(chr, start, chr, end, rTreeOffset);

            Integer chrIdx = reader.getIdForChr(chr);

            List<LocusScore> features = new ArrayList<>();
            for (byte[] c : chunks) {
                reader.decodeWigData(c, chrIdx, start, end, features);
            }

            final int size = features.size();
            int[] starts = new int[size];
            int[] ends = new int[size];
            float[] values = new float[size];
            for (int i = 0; i < size; i++) {
                final LocusScore locusScore = features.get(i);
                starts[i] = locusScore.getStart();
                ends[i] = locusScore.getEnd();
                values[i] = locusScore.getScore();
            }
            return new DataTile(starts, ends, values, null);

        } catch (IOException e) {
            throw new RuntimeException(e);
        }

    }


    /**
     * Return bigwig "zoom data" if available for the resolution encoded by "zoom"
     *
     * @param chr
     * @param start
     * @param end
     * @param zoom
     * @return
     */
    @Override
    protected List<LocusScore> getPrecomputedSummaryScores(String chr, int start, int end, int zoom) {
        // Translate IGV zoom number to bpperpixel.
        try {
            if (Globals.CHR_ALL.equals(chr)) {
                if (genome.getHomeChromosome().equals(Globals.CHR_ALL) && windowFunction != WindowFunction.none) {
                    return this.getWholeGenomeScores();
                } else {
                    return null;
                }
            }


            Chromosome chromosome = genome.getChromosome(chr);
            if (chromosome == null) {
                throw new RuntimeException("Unexpected chromosome name: " + chr);
            }


            double nBins = Math.pow(2, zoom);
            double scale = chromosome.getLength() / (nBins * 700);
            BBZoomHeader zlHeader = reader.zoomLevelForScale(scale);

            if (zlHeader == null) {
                return null;
            } else {
                long rTreeOffset = zlHeader.indexOffset;
                int chrIdx = reader.getIdForChr(chr);
                List<byte[]> chunks = this.reader.getLeafChunks(chr, start, chr, end, rTreeOffset);
                List<LocusScore> features = new ArrayList<>();
                for (byte[] c : chunks) {
                    reader.decodeZoomData(c, chrIdx, start, end, windowFunction, features);
                }
                return features;
            }
        } catch (IOException e) {
            // Wrap the IOException for now, to avoid refactoring of hierarchy
            throw new RuntimeException(e);
        }
    }

    @Override
    public int getLongestFeature(String chr) {
        return 0;
    }


    public TrackType getTrackType() {
        return TrackType.OTHER;
    }


    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return availableWindowFunctions;
    }

    @Override
    public void dispose() {

    }


    private List<LocusScore> getWholeGenomeScores() {

        try {
            if (genome.getHomeChromosome().equals(Globals.CHR_ALL) && windowFunction != WindowFunction.none) {

                if (!wholeGenomeScores.containsKey(windowFunction)) {

                    int screenWidth = 1000;  // nominal
                    double scale = genome.getWGLength() / screenWidth;

                    int maxChromId = reader.getChromosomeNames().length - 1;
                    String firstChr = reader.getChrForId(0);
                    String lastChr = reader.getChrForId(maxChromId);

                    ArrayList<LocusScore> scores = new ArrayList<LocusScore>();
                    wholeGenomeScores.put(windowFunction, scores);

                    BBZoomHeader lowestResHeader = reader.zoomLevelForScale(scale, 1000);
                    if (lowestResHeader == null) return null;

                    Set<String> wgChrNames = new HashSet<>(genome.getLongChromosomeNames());

                    long rTreeOffset = lowestResHeader.indexOffset;
                    List<byte[]> chunks = this.reader.getLeafChunks(firstChr, 0, lastChr, Integer.MAX_VALUE, rTreeOffset);

                    List<LocusScore> features = new ArrayList<>();
                    for (byte[] c : chunks) {
                        reader.decodeZoomData(c, -1, -1, -1, windowFunction, features);
                    }

                    for (LocusScore s : features) {
                        WigDatum rec = (WigDatum) s;   // Unfortunate, but this is java
                        String chr = genome.getCanonicalChrName(rec.getChr());
                        if (wgChrNames.contains(chr)) {
                            int genomeStart = genome.getGenomeCoordinate(chr, rec.getStart());
                            int genomeEnd = genome.getGenomeCoordinate(chr, rec.getEnd());
                            scores.add(new BasicScore(genomeStart, genomeEnd, rec.getScore()));
                        }
                    }
                    scores.sort(Comparator.comparingInt(Locatable::getStart));

                }
                return wholeGenomeScores.get(windowFunction);
            } else {
                return null;
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

}
