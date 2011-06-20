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
import org.broad.igv.data.AbstractDataSource;
import org.broad.igv.data.BasicScore;
import org.broad.igv.data.DataTile;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.util.collections.FloatArrayList;
import org.broad.igv.util.collections.IntArrayList;
import org.broad.tribble.util.SeekableStream;
import org.broad.tribble.util.SeekableStreamFactory;

import java.io.IOException;
import java.util.*;

/**
 * @author jrobinso
 * @date Jun 19, 2011
 */
public class BigWigDataSource extends AbstractDataSource {

    Collection<WindowFunction> availableWindowFunctions =
            Arrays.asList(WindowFunction.min, WindowFunction.mean, WindowFunction.max);
    WindowFunction windowFunction = WindowFunction.mean;

    BBFileReader reader;
    private BBZoomLevels levels;
    String path;
    SeekableStream ss;

    public BigWigDataSource(String path, Genome genome) throws IOException {
        super(genome);
        this.path = path;

        ss = SeekableStreamFactory.getStreamFor(path);
        reader = new BBFileReader(path, ss);

        if (reader.isBigBedFile()) {
            throw new RuntimeException("BigBed files are not currently supported (coming soon)");
        }

        levels = reader.getZoomLevels();
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
        Chromosome c = genome.getChromosome(chr);
        if (c == null) {
            return null;
        }
        int l = c.getLength();
        double scale =  l / (Math.pow(2, zoom) * 700);   // TODO -- use actual screen resolution

        // Find first zoom level

        //
        int bbLevel = -1;
        int reductionLevel = -1;
        for (BBZoomLevelHeader zlHeader : levels.getZoomLevelHeaders()) {
            reductionLevel = zlHeader.getReductionLevel();
            bbLevel = zlHeader.getZoomLevel();
            if (reductionLevel > scale) {
                break;
            }
        }

        // If we are at the highest precomputed resolution compare to the requested resolution.  If they differ
        // by more than a factor of 2 compute "on the fly"

        if (bbLevel >1 || (bbLevel == 1 && (reductionLevel / scale) < 2)) {
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
