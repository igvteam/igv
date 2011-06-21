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

package org.broad.igv.bigwig;

import org.broad.igv.Globals;
import org.broad.igv.bbfile.*;
import org.broad.igv.data.*;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.FeatureSource;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.util.collections.FloatArrayList;
import org.broad.igv.util.collections.IntArrayList;
import org.broad.tribble.util.SeekableStream;
import org.broad.tribble.util.SeekableStreamFactory;

import java.io.IOException;
import java.util.*;

/**
 * A hybrid source, implements both DataSource and FeatureSource.
 *
 * @author jrobinso
 * @date Jun 19, 2011
 */
public class BigWigDataSource extends AbstractDataSource implements FeatureSource {

    Collection<WindowFunction> availableWindowFunctions =
            Arrays.asList(WindowFunction.min, WindowFunction.mean, WindowFunction.max);
    WindowFunction windowFunction = WindowFunction.mean;

    BBFileReader reader;
    private BBZoomLevels levels;

    // Feature visibility window (for bigBed)
    int featureVisiblityWindow = -1;
    private GenomeSummaryData genomeSummaryData;
    private List<LocusScore> wholeGenomeScores;


    public BigWigDataSource(String path, Genome genome) throws IOException {
        super(genome);

        SeekableStream ss = SeekableStreamFactory.getStreamFor(path);
        reader = new BBFileReader(path, ss);
        levels = reader.getZoomLevels();
    }

    public BigWigDataSource(BBFileReader reader, Genome genome) throws IOException {
        super(genome);

        this.reader = reader;
        levels = reader.getZoomLevels();

        // Assume 1000 pixel screen, pick visibility level to be @ highest resolution zoom.
        // TODO -- something smarter, like scaling by actual density
        if (levels.getZoomHeaderCount() > 0) {
            BBZoomLevelHeader firstLevel = levels.getZoomLevelHeaders().get(0);
            featureVisiblityWindow = firstLevel.getReductionLevel() * 2000;
        }
    }


    public double getDataMax() {
        return 100;
    }

    public double getDataMin() {
        return 0;
    }

    public TrackType getTrackType() {
        return TrackType.OTHER;
    }

    public void setWindowFunction(WindowFunction statType) {
        this.windowFunction = statType;
    }

    public boolean isLogNormalized() {
        return false;
    }

    public void refreshData(long timestamp) {

    }

    @Override
    public int getLongestFeature(String chr) {
        return 0;
    }

    public WindowFunction getWindowFunction() {
        return windowFunction;
    }

    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return availableWindowFunctions;
    }

    @Override
    protected List<LocusScore> getPrecomputedSummaryScores(String chr, int start, int end, int zoom) {

        if (chr.equals(Globals.CHR_ALL)) {
            return getWholeGenomeScores();
        } else {
            return getZoomSummaryScores(chr, start, end, zoom);
        }
    }


    private BBZoomLevelHeader getZoomLevelForScale(double scale) {
        final ArrayList<BBZoomLevelHeader> headers = levels.getZoomLevelHeaders();

        for (BBZoomLevelHeader zlHeader : headers) {
            int reductionLevel = zlHeader.getReductionLevel();
            if (reductionLevel > scale) {
                return zlHeader;
            }
        }
        return headers.get(headers.size() - 1);

    }

    protected List<LocusScore> getZoomSummaryScores(String chr, int start, int end, int zoom) {

        Chromosome c = genome.getChromosome(chr);
        if (c == null) {
            return null;
        }
        int l = c.getLength();
        double scale = l / (Math.pow(2, zoom) * 700);   // TODO -- use actual screen resolution


        //
        BBZoomLevelHeader zlHeader = getZoomLevelForScale(scale);
        int bbLevel = zlHeader.getZoomLevel();
        int reductionLevel = zlHeader.getReductionLevel();


        // If we are at the highest precomputed resolution compare to the requested resolution.  If they differ
        // by more than a factor of 2 compute "on the fly"

        if (reader.isBigBedFile() || bbLevel > 1 || (bbLevel == 1 && (reductionLevel / scale) < 2)) {
            ArrayList<LocusScore> scores = new ArrayList(1000);
            ZoomLevelIterator zlIter = reader.getZoomLevelIterator(bbLevel, chr, start, chr, end, false);
            while (zlIter.hasNext()) {
                ZoomDataRecord rec = zlIter.next();

                float v;
                switch (windowFunction) {
                    case min:
                        v = rec.getMinVal();
                        break;
                    case max:
                        v = rec.getMaxVal();
                        break;
                    default:
                        v = rec.getMeanVal();

                }

                BasicScore bs = new BasicScore(rec.getChromStart(), rec.getChromEnd(), v);
                scores.add(bs);
            }
            return scores;

        } else {
            // No precomputed scores for this resolution level
            return null;
        }
    }


    RawDataInterval currentInterval = null;


    @Override
    protected synchronized DataTile getRawData(String chr, int start, int end) {

        if (chr.equals(Globals.CHR_ALL)) {
            return null;
        }


        if (currentInterval != null && currentInterval.contains(chr, start, end)) {
            return currentInterval.tile;
        }

        // TODO -- fetch data directly in arrays to avoid creation of multiple "WigItem" objects?
        IntArrayList startsList = new IntArrayList(100000);
        IntArrayList endsList = new IntArrayList(100000);
        FloatArrayList valuesList = new FloatArrayList(100000);

        Iterator<WigItem> iter = reader.getBigWigIterator(chr, start, chr, end, false);

        while (iter.hasNext()) {
            WigItem wi = iter.next();
            startsList.add(wi.getStartBase());
            endsList.add(wi.getEndBase());
            valuesList.add(wi.getWigValue());
        }

        DataTile tile = new DataTile(startsList.toArray(), endsList.toArray(), valuesList.toArray(), null);
        currentInterval = new RawDataInterval(chr, start, end, tile);

        return tile;

    }

    private List<LocusScore> getWholeGenomeScores() {

        if (genome.getHomeChromosome().equals(Globals.CHR_ALL)) {
            if (wholeGenomeScores == null) {
                wholeGenomeScores = new ArrayList<LocusScore>();
                for (Chromosome chr : genome.getChromosomes()) {

                    double scale = chr.getLength() / 1000;  // TODO - use screen width instead of constant "1000"
                    BBZoomLevelHeader lowestResHeader = this.getZoomLevelForScale(scale);

                    int lastGenomeEnd = -1;
                    String chrName = chr.getName();
                    int end = chr.getLength();

                    ZoomLevelIterator zlIter = reader.getZoomLevelIterator(
                            lowestResHeader.getZoomLevel(), chrName, 0, chrName, end, false);
                    while (zlIter.hasNext()) {
                        ZoomDataRecord rec = zlIter.next();
                        int genomeStart = genome.getGenomeCoordinate(chrName, rec.getChromStart());
                        if (genomeStart < lastGenomeEnd) {
                            continue;
                        }

                        int genomeEnd = genome.getGenomeCoordinate(chrName, rec.getChromEnd());
                        float value = rec.getMeanVal();  // TODO window function
                        wholeGenomeScores.add(new BasicScore(genomeStart, genomeEnd, value));
                        lastGenomeEnd = genomeEnd;
                    }
                }

            }
            return wholeGenomeScores;
        } else {
            return null;
        }

    }


    // Feature interface follows ------------------------------------------------------------------------

    public Iterator getFeatures(String chr, int start, int end) throws IOException {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {
        return this.getSummaryScoresForRange(chr, start, end, zoom);
    }

    public int getFeatureWindowSize() {
        return this.featureVisiblityWindow;
    }

    public void setFeatureWindowSize(int size) {
        this.featureVisiblityWindow = size;
    }

    public Class getFeatureClass() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }


    //  End FeatureSource interface ----------------------------------------------------------------------

    static class RawDataInterval {
        String chr;
        int start;
        int end;
        DataTile tile;

        RawDataInterval(String chr, int start, int end, DataTile tile) {
            this.chr = chr;
            this.start = start;
            this.end = end;
            this.tile = tile;
        }

        public boolean contains(String chr, int start, int end) {
            return chr.equals(this.chr) && start >= this.start && end <= this.end;
        }
    }


}
